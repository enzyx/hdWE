#!/usr/bin/python2
"""
hdWE is a hyperdimensional weighted ensemble simulation implementation.
"""
from __future__ import print_function
import sys
import ConfigParser
import argparse
import numpy
import json
import initiate
from iteration import Iteration
from logger import Logger
import threading
import reweighting
from thread_container import ThreadContainer
from hdWE_parameters import HdWEParameters
from analysis_operations import meanRateMatrix

###### Parse command line ###### 

parser =argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-c', '--conf', type=str, dest="configfile", 
                    metavar="FILE", required=False, default='hdWE.conf',
                    help="The hdWE configuration file")      
parser.add_argument("--append", dest="append", action='store_true', 
                    default=False,
                    help="continue previous iterations with parameters from conf file.")  
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
hdWE_parameters.loadConfParameters(config, args.debug)

# Initialize the logger
logger = Logger(filename=hdWE_parameters.logfile, 
                append=args.append,
                debug=args.debug)
logger.logParameters(hdWE_parameters)
if args.append:
    iterations = logger.loadIterations()  

# Setup the workdir and initiate iterations
initiate.prepare(hdWE_parameters.workdir, hdWE_parameters.starting_structure, args.append, args.debug)
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
    iteration = Iteration(iteration_counter)
    parent_iteration = iterations[iteration_counter - 1]
    # Generate all previous bins for new iteration
    for parent_bin in parent_iteration:
        iteration.generateBin(reference_iteration_id=parent_bin.getReferenceIterationId(),
                              reference_bin_id=parent_bin.getReferenceBinId(),
                              reference_segment_id=parent_bin.getReferenceSegmentId(),
                              target_number_of_segments=hdWE_parameters.segments_per_bin)
    
    coordinates = numpy.array([]) # numpy array?
    
    max_bin_probability = parent_iteration.getMaxBinProbability()
    
    # Sort segments in bins. Generate new bin if required
    for parent_bin in parent_iteration:
        for segment in parent_bin:
            coordinates = md_module.calcRmsdToBins(segment, iteration.bins)
            min_coordinate = numpy.min(coordinates)
            # Sort Segment into appropriate bin
            if (min_coordinate <= hdWE_parameters.coordinate_threshold or 
               segment.getProbability() < hdWE_parameters.minimal_probability*max_bin_probability or 
               iteration.getNumberOfBins() >= hdWE_parameters.max_bins):
                   
                bin_id = numpy.argmin(coordinates)
                iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                  parent_bin_id=segment.getBinId(),
                                                  parent_segment_id=segment.getId())
            # If necessary create new bin
            else:
                bin_id = iteration.generateBin(reference_iteration_id=segment.getIterationId(),
                                      reference_bin_id=segment.getBinId(),
                                      reference_segment_id=segment.getId(),
                                      target_number_of_segments=hdWE_parameters.segments_per_bin)
                iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                  parent_bin_id=segment.getBinId(),
                                                  parent_segment_id=segment.getId())
    # Split and merge (Manu)
    # Parallel
    if config.get('hdWE','number-of-threads') > 1:
        thread_container = ThreadContainer()
        for this_bin in iteration:
            thread_container.appendJob(threading.Thread(target=this_bin.resampleSegments))
            if thread_container.getNumberOfJobs() >= config.get('hdWE','number-of-threads'):
                thread_container.runJobs()
        # Run remaining jobs
        thread_container.runJobs()
    # Serial
    else:
        for this_bin in iteration:
            this_bin.resampleSegments()

    if args.debug: 
        print("The overall probabiliy is {0:05f}".format(iteration.getProbability()))
        
    # check which bins have to be rerun
    reweighting.checkOutratesForConvergence(iterations, hdWE_parameters.reweighting_range, 0.05)
                        
				

    # Run MD
    md_module.RunMDs(iteration)

    iterations.append(iteration)

    #reweight Bin Probabilities
    if iteration.getNumberOfBins() > 1:
        if hdWE_parameters.reweighting_range > 0:
            reweighting.reweightBinProbabilities(iterations, hdWE_parameters.reweighting_range)

    
    # log iteration (Rainer)
    logger.logIteration(iteration)
    
logger.close()

#count total n of segments
n_segments = 0
for iteration_loop in iterations:
     n_segments += iteration_loop.getNumberOfSegments()
    
print('hdWE completed. Total number of propagated segments: ' + str(n_segments) + '            ')
