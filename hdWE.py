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
import os
import constants

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

#############################
#        Variables          #
#############################
CONFIGFILE = args.configfile
config = ConfigParser.ConfigParser()
config.read(CONFIGFILE)

WORKDIR               = config.get('hdWE','WORKDIR')
JOBNAME               = config.get('hdWE','JOBNAME')
LOGDIR                = constants.getLogDirPath(WORKDIR, JOBNAME)
RUNDIR                = constants.getRunDirPath(WORKDIR, JOBNAME)
APPEND                = args.append
OVERWRITE             = args.overwrite
DEBUG                 = args.debug
MAX_ITERATIONS        = int(config.get('hdWE','max-iterations'))
SEGMENTS_PER_BIN      = int(config.get('hdWE','segments-per-bin'))
STARTING_STRUCTURE    = config.get('hdWE','starting-structure')
COORDINATE_THRESHOLD  = float(config.get('hdWE','threshold'))
MINIMAL_PROBABILITY   = float(config.get('hdWE','minimal-probability'))
MAX_BINS              = int(config.get('hdWE','max-bins'))
REWEIGHTING_RANGE     = float(config.get('hdWE','reweighting-range'))
CONVERGENCE_RANGE     = float(config.get('hdWE','convergence-range'))
CONVERGENCE_THRESHOLD = float(config.get('hdWE','convergence-threshold'))
NUMBER_OF_THREADS     = int(config.get('hdWE','number-of-threads'))
MERGE_MODE            = int(config.get('hdWE','merge-mode'))

if "amber" in config.sections():
    MD_PACKAGE = "amber"
elif "gromacs" in config.sections():
    MD_PACKAGE = "gromacs"
else:
    raise Exception("No MD package (amber/gromacs) section in configuration file")

# The global list of arrays
iterations = [] 

#############################
#        Check files        #
#############################
if not os.path.isfile(CONFIGFILE):
    print("Configuration file not found: " + CONFIGFILE)
    sys.exit(-1)

#############################
#           Main            #
#############################
# Setup the WORKDIR
initiate.prepare(WORKDIR,
                 JOBNAME, 
                 STARTING_STRUCTURE, 
                 OVERWRITE, 
                 APPEND, 
                 DEBUG)

# Initialize the logger
logger = Logger(LOGDIR)

# Initiate iterations
if APPEND:
    iterations = logger.loadIterations()
    #TODO: check if all files are present  
else:
    iterations.append(initiate.create_initial_iteration(SEGMENTS_PER_BIN))
    logger.log(iterations[0], CONFIGFILE)

# Check MD suite
if(MD_PACKAGE == "amber"):
    from amber_module import MD_module
    md_module = MD_module(CONFIGFILE, DEBUG)
if(MD_PACKAGE == "gromacs"):
    print("Sorry, support for gromacs is not implemented yet.")
    sys.exit(-1)


#############################
#         Main Loop         #
#############################
for iteration_counter in range(len(iterations), MAX_ITERATIONS + 1):
    sys.stdout.write('\033[1m\n hdWE Status:\033[0m ' + 'Iteration ' + str(iteration_counter).zfill(5) + '\n')    
    sys.stdout.flush()

    # 1. Initialization of new iteration. 
    #    - initialize iteration and parent_iteration.
    #    - copy existing bins from the previous iteration into the new iteration.
    sys.stdout.write(' Initialization\n')
    sys.stdout.flush()

    iteration = Iteration(iteration_counter)
    parent_iteration = iterations[iteration_counter - 1]
    for parent_bin in parent_iteration:
        iteration.generateBin(reference_iteration_id=parent_bin.getReferenceIterationId(),
                              reference_bin_id=parent_bin.getReferenceBinId(),
                              reference_segment_id=parent_bin.getReferenceSegmentId(),
                              target_number_of_segments=SEGMENTS_PER_BIN,
                              outrates_converged = parent_bin.isConverged())

    # 2. Sorting of segments into bins
    #    - Calculate the rmsd of each segment to all existing bins.
    #    - Generate new bin if required
    sys.stdout.write(' Sorting Segments into Bins\n')
    sys.stdout.flush()    
    coordinates = numpy.array([])
    max_bin_probability = parent_iteration.getMaxBinProbability()
    rmsd_matrix = md_module.calcRmsdSegmentsToBinsMatrix(parent_iteration)

    segment_id = 0
    new_bins   = []
    for parent_bin in parent_iteration:
        for segment in parent_bin:
            is_segment_handled = False
            min_coordinate = numpy.min(rmsd_matrix[segment_id])
            min_target_bin_id = numpy.argmin(rmsd_matrix[segment_id])
            first_coordinate_below_threshold = COORDINATE_THRESHOLD + 1.0
            first_target_bin_id = -1
            # find first bin below threshold:
            for target_bin_id, coordinate in enumerate(rmsd_matrix[segment_id]):
                if coordinate <= COORDINATE_THRESHOLD:
                    first_coordinate_below_threshold = coordinate
                    first_target_bin_id = target_bin_id
                    break
            min_coordinate = first_coordinate_below_threshold
            min_target_bin_id = first_target_bin_id
            # Sort Segment into appropriate bin
            if (min_coordinate <= COORDINATE_THRESHOLD or 
               segment.getProbability() < MINIMAL_PROBABILITY*max_bin_probability or 
               iteration.getNumberOfBins() >= MAX_BINS):
                bin_id = min_target_bin_id
                iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                  parent_bin_id=segment.getBinId(),
                                                  parent_segment_id=segment.getId())
                is_segment_handled = True
            # If necessary create new bin
            # Check if this segment fits into one of the new bins first
            if not is_segment_handled and len(new_bins) > 0:
                new_rmsds = md_module.calcRmsdToBins(segment, new_bins)
                min_coordinate = numpy.min(new_rmsds)
                first_coordinate_below_threshold = COORDINATE_THRESHOLD + 1.0
                first_target_bin_id = -1
                # find first bin below threshold:
                for target_bin_id, coordinate in enumerate(new_rmsds):
                    if coordinate <= COORDINATE_THRESHOLD:
                        first_coordinate_below_threshold = coordinate
                        first_target_bin_id = target_bin_id
                        break
                min_coordinate = first_coordinate_below_threshold
                min_target_bin_id = first_target_bin_id
                if min_coordinate <= COORDINATE_THRESHOLD:
                    bin_id = min_target_bin_id
                    iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                parent_bin_id=segment.getBinId(),
                                                parent_segment_id=segment.getId())
                    is_segment_handled = True
                    if DEBUG:
                        print("Segment {0} fits in new bin {1}".format(
                                    segment.getNameString(), 
                                    new_bins[min_target_bin_id]))
            if not is_segment_handled:
                # Ok seems that we need a new bin
                bin_id = iteration.generateBin(reference_iteration_id=segment.getIterationId(),
                                    reference_bin_id=segment.getBinId(),
                                    reference_segment_id=segment.getId(),
                                    target_number_of_segments=SEGMENTS_PER_BIN)
                iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                parent_bin_id=segment.getBinId(),
                                                parent_segment_id=segment.getId())
                new_bins.append(iteration.bins[bin_id])
                is_segment_handled = True
                if DEBUG:
                    print("created new bin {0} from segment ({1})".format(
                                    new_bins[-1].getId(), 
                                    segment.getNameString()))
            segment_id += 1
    
    # 3. Backup the segments lists of all bins
    #    Saving the segment assignments for correct
    #    rate matrix calculation in the original_segments
    #    list of bins. This list should be immutable 
    for _bin in iteration:
        _bin.backupInitialSegments()
    
    #Append iteration that has now all bins here!        
    iterations.append(iteration)            
            
    # 4. Reweighting of bin probabilities
    #    The order of the following steps should no longer matter.  
    if iteration.getNumberOfBins() > 1 and REWEIGHTING_RANGE > 0.0:
        sys.stdout.write(' Reweighting Bin Probabilities\n')
        sys.stdout.flush()
        #TODO: why does reweighting need the workdir and jobname? fix it! 
        reweighting.reweightBinProbabilities(iterations,
                                             REWEIGHTING_RANGE,
                                             WORKDIR,
                                             JOBNAME)
    # 5. Check convergence of the bins
    if iteration_counter > CONVERGENCE_RANGE:    
        sys.stdout.write(' Check Bin Convergence\n')
        sys.stdout.flush() 
        convergenceCheck.checkOutratesForConvergence(iterations, 
                                                     CONVERGENCE_RANGE,
                                                     CONVERGENCE_THRESHOLD,
                                                     DEBUG)

    # 5. Resampling
    sys.stdout.write(' Resampling\n')
    sys.stdout.flush() 
    # Parallel
    if NUMBER_OF_THREADS > 1:
        thread_container = ThreadContainer()
        for this_bin in iterations[-1]:
            thread_container.appendJob(threading.Thread(target=this_bin.resampleSegments(MERGE_MODE) ))
            if thread_container.getNumberOfJobs() >= NUMBER_OF_THREADS:
                thread_container.runJobs()
        # Run remaining jobs
        thread_container.runJobs()
    # Serial
    else:
        for this_bin in iterations[-1]:
            this_bin.resampleSegments(MERGE_MODE)

    # 6. Run MDs
    sys.stdout.write(' Run MDs\n')
    sys.stdout.flush() 
    md_module.RunMDs(iterations[-1])
    
    
    # 7. log everything
    logger.log(iterations[-1], CONFIGFILE)

    #if DEBUG: 
    print(" The overall probabiliy is {0:05f}".format(iterations[-1].getProbability()))

    #count total n of segments during iterations
    n_segments = 0
    for iteration_loop in iterations:
        n_segments += iteration_loop.getNumberOfPropagatedSegments()
    sys.stdout.write(' (Total number of propagated segments: ' + str(n_segments-1)+ ')\n')
    sys.stdout.flush()

    all_converged = True
    for bin_loop in iterations[-1]:
        if bin_loop.isConverged() == False:
            all_converged = False
    if all_converged == True:
        print('\nExploring Mode: All bins are converged                                         ')
        break        
    #check for empty bins
    empty_bins = 0
    for bin_loop in iterations[-1]:
        if bin_loop.getNumberOfSegments() == 0:
            empty_bins += 1
    print('    Empty bins: ' + str(empty_bins))

#############################
#     End of Main Loop      #
#############################

#count total n of segments
n_segments = 0
for iteration_loop in iterations:
    n_segments += iteration_loop.getNumberOfPropagatedSegments()
print('\nhdWE completed. Total number of propagated segments: ' + str(n_segments-1) + '            ')



