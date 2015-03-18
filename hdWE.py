#!/usr/bin/python2
"""
hdWE is a hyperdimensional weighted ensemble simulation implementation.
"""
from __future__ import print_function
import sys
import ConfigParser
import argparse
import numpy
import initiate
from iteration import Iteration
from logger import Logger
import threading
import reweighting
import convergenceCheck
from thread_container import ThreadContainer
from hdWE_parameters import HdWEParameters

###### Parse command line ###### 

parser =argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-c', '--conf', type=str, dest="configfile", 
                    metavar="FILE", required=False, default='hdWE.conf',
                    help="The hdWE configuration file")      
prev_data = parser.add_mutually_exclusive_group()
prev_data.add_argument("--append", dest="append", action='store_true', 
                    default=False,
                    help="continue previous iterations with parameters from conf file.")  
prev_data.add_argument("--overwrite", dest="overwrite", action='store_true', 
                    default=False,
                    help="overwrite previous simulation data.")  
parser.add_argument('--debug', dest="debug", action="store_true",
                    default=False, help="Turn debugging on")

args = parser.parse_args()                

###### Parse Config File #######
config = ConfigParser.ConfigParser()
config.read(args.configfile)

#############################
# Main
#############################

# The global list of arrays
iterations = [] 

# Initialize the parameter container
hdWE_parameters = HdWEParameters()
hdWE_parameters.loadConfParameters(config, args.configfile, args.debug)

# Setup the workdir
initiate.prepare(hdWE_parameters.workdir,
                 hdWE_parameters.jobname, 
                 hdWE_parameters.starting_structure, 
                 args.overwrite, 
                 args.append, 
                 args.debug)


# Initialize the logger
logger = Logger(filename="{jn}.log".format(jn=hdWE_parameters.jobname), 
                append=args.append,
                debug=args.debug)
logger.logParameters(hdWE_parameters)
if args.append:
    iterations = logger.loadIterations()  

# initiate iterations
if len(iterations)==0:
    iterations.append(initiate.create_initial_iteration(hdWE_parameters.segments_per_bin))
    logger.logIteration(iterations[0])


# Check MD suite
if(hdWE_parameters.md_package.lower() == "amber"):
    from amber_module import MD_module
    md_module = MD_module(args.configfile, debug=args.debug)
if(hdWE_parameters.md_package.lower() == "gromacs"):
    print("Sorry, support for gromacs is not implemented yet.")
    sys.exit(-1)


# Loop 
for iteration_counter in range(len(iterations), hdWE_parameters.max_iterations + 1):
    # 1. Initialisation of new iteration. 
    #    - initialize iteration and parent_iteration.
    #    - copy existing bins from the previous iteration into the new iteration.
    iteration = Iteration(iteration_counter)
    parent_iteration = iterations[iteration_counter - 1]
    for parent_bin in parent_iteration:
        iteration.generateBin(reference_iteration_id=parent_bin.getReferenceIterationId(),
                              reference_bin_id=parent_bin.getReferenceBinId(),
                              reference_segment_id=parent_bin.getReferenceSegmentId(),
                              target_number_of_segments=hdWE_parameters.segments_per_bin,
                              outrates_converged = parent_bin.isConverged())

    # 2. Sorting of segments into bins
    #    - Calculate the rmsd of each segment to all existing bins.
    #    - Generate new bin if required

    coordinates = numpy.array([])
    max_bin_probability = parent_iteration.getMaxBinProbability()
    rmsd_matrix = md_module.calcRmsdSegmentsToBinsMatrix(parent_iteration)
    segment_id = 0
    new_bins   = []
    for parent_bin in parent_iteration:
        for segment in parent_bin:
            min_coordinate = numpy.min(rmsd_matrix[segment_id])
            # Sort Segment into appropriate bin
            if (min_coordinate <= hdWE_parameters.coordinate_threshold or 
               segment.getProbability() < hdWE_parameters.minimal_probability*max_bin_probability or 
               iteration.getNumberOfBins() >= hdWE_parameters.max_bins):
                bin_id = numpy.argmin(rmsd_matrix[segment_id])
                iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                  parent_bin_id=segment.getBinId(),
                                                  parent_segment_id=segment.getId())
            # If necessary create new bin
            else:
                # Check if this segment fits into one of the new bins first
                is_segment_handled = False
                if len(new_bins) > 0:
                    rmsds = md_module.calcRmsdToBins(segment, new_bins)
                    min_coordinate = numpy.min(rmsds)
                    if min_coordinate <= hdWE_parameters.coordinate_threshold:
                        bin_id = new_bins[numpy.argmin(rmsds)].getId()
                        iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                  parent_bin_id=segment.getBinId(),
                                                  parent_segment_id=segment.getId())
                        is_segment_handled = True
                        if args.debug:
                            print("Segment {0} fits in new bin {1}".format(
                                        segment.getNameString(), 
                                        new_bins[numpy.argmin(rmsds)].getId()))
                if not is_segment_handled:
                    # Ok seems that we need a new bin
                    bin_id = iteration.generateBin(reference_iteration_id=segment.getIterationId(),
                                      reference_bin_id=segment.getBinId(),
                                      reference_segment_id=segment.getId(),
                                      target_number_of_segments=hdWE_parameters.segments_per_bin)
                    iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                  parent_bin_id=segment.getBinId(),
                                                  parent_segment_id=segment.getId())
                    new_bins.append(iteration.bins[bin_id])
                    if args.debug:
                        print("created new bin {0} from segment ({1})".format(
                                        new_bins[-1].getId(), 
                                        segment.getNameString()))
            segment_id += 1
    
    #Append iteration that has now all bins here!        
    iterations.append(iteration)            
            
    # 3. Reweighting of bin probabilities
    #    (Has to happen before resampling)
    if iteration.getNumberOfBins() > 3 and hdWE_parameters.reweighting_range > 0.0:
            reweighting.reweightBinProbabilities(iterations,
                                                 hdWE_parameters.reweighting_range,
                                                 hdWE_parameters.workdir,
                                                 hdWE_parameters.jobname)
    # 4. Check convergence of the bins
    if iteration_counter > hdWE_parameters.convergence_range:    
        convergenceCheck.checkOutratesForConvergence(iterations, 
                                                     hdWE_parameters.convergence_range,
                                                     hdWE_parameters.convergence_threshold,
                                                     args.debug)




    # 5. Resampling
    # Parallel
    if config.get('hdWE','number-of-threads') > 1:
        thread_container = ThreadContainer()
        for this_bin in iterations[-1]:
            thread_container.appendJob(threading.Thread(target=this_bin.resampleSegments))
            if thread_container.getNumberOfJobs() >= config.get('hdWE','number-of-threads'):
                thread_container.runJobs()
        # Run remaining jobs
        thread_container.runJobs()
    # Serial
    else:
        for this_bin in iterations[-1]:
            this_bin.resampleSegments()

    # 6. Run MDs
    md_module.RunMDs(iterations[-1])
    
    
    # 7. log iteration
    logger.logIteration(iterations[-1])


    #if args.debug: 
    print("The overall probabiliy is {0:05f}".format(iterations[-1].getProbability()))

    #count total n of segments during iterations
    n_segments = 0
    for iteration_loop in iterations:
        n_segments += iteration_loop.getNumberOfPropagatedSegments()
    sys.stdout.write(' (Total number of propagated segments: ' + str(n_segments)+ ')')
    sys.stdout.flush()

    all_converged = True
    for bin_loop in iterations[-1]:
        if bin_loop.isConverged() == False:
            all_converged = False
    if all_converged == True:
        print('\nExploring Mode: All bins are converged                                         ')
        break        

# END of iteration loop

#count total n of segments
n_segments = 0
for iteration_loop in iterations:
    n_segments += iteration_loop.getNumberOfPropagatedSegments()
print('\nhdWE completed. Total number of propagated segments: ' + str(n_segments) + '            ')
   
logger.close()


