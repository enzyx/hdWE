#!/usr/bin/python2
"""
hdWE is a hyperdimensional weighted ensemble simulation implementation.
"""
from __future__ import print_function
import sys
import ConfigParser
import argparse
import initiate
from iteration import Iteration
from logger import Logger
import threading
import resorting
import recycling
import reweighting
from thread_container import ThreadContainer
import os
import glob
import constants

#### Parse command line #### 

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-c', '--conf', type=str, dest="configfile", 
                    metavar="FILE", required=False, default='hdWE.conf',
                    help="The hdWE configuration file")
parser.add_argument('-d', '--debug', dest="debug", action="store_true",
                    default=False, help="Turn debugging on")
prev_data = parser.add_mutually_exclusive_group()
prev_data.add_argument('-a', '--append', dest="append", action='store_true', 
                       default=False,
                       help='continue previous iterations with parameters from conf file')  
prev_data.add_argument('-o', '--overwrite', dest="overwrite", action='store_true', 
                       default=False,
                       help='overwrite previous simulation data')
                    
args = parser.parse_args()                

#############################
#        Variables          #
#############################
CONFIGFILE = args.configfile
config = ConfigParser.ConfigParser()
config.read(CONFIGFILE)

WORKDIR                           = config.get('hdWE','WORKDIR')
JOBNAME                           = config.get('hdWE','JOBNAME')
LOGDIR                            = constants.getLogDirPath(WORKDIR, JOBNAME)
RUNDIR                            = constants.getRunDirPath(WORKDIR, JOBNAME)
APPEND                            = args.append
OVERWRITE                         = args.overwrite
DEBUG                             = args.debug
MAX_ITERATIONS                    = int(config.get('hdWE','max-iterations'))
INITIAL_TARGET_NUMBER_OF_SEGMENTS = int(config.get('hdWE','segments-per-bin'))
INITIAL_BOUNDARIES                = initiate.parseInitialBoundaries(config)
STARTING_STRUCTURE                = config.get('hdWE','starting-structure')
REWEIGHTING_RANGE                 = float(config.get('hdWE','reweighting-range'))
NUMBER_OF_THREADS                 = int(config.get('hdWE','number-of-threads'))
MERGE_MODE                        = int(config.get('hdWE','merge-mode'))
STEADY_STATE                      = bool(config.get('hdWE', 'steady-state').lower() == "true")
START_STATES                      = initiate.parseState(config.get('hdWE', 'start-state'))
END_STATES                        = initiate.parseState(config.get('hdWE', 'end-state'))

if "amber" in config.sections():
    MD_PACKAGE = "amber"
elif "gromacs" in config.sections():
    MD_PACKAGE = "gromacs"
else:
    raise Exception("No MD package (amber/gromacs) section in configuration file")

# The global list of iterations
iterations = [] 

#############################
#        Check files        #
#############################
if not os.path.isfile(CONFIGFILE):
    print("Configuration file not found: " + CONFIGFILE)
    sys.exit(-1)
if APPEND and not (os.path.exists(LOGDIR) and
                   os.path.exists(RUNDIR) and 
                   glob.glob(LOGDIR + "*.iter")):
    sys.stderr.write("WARNING: appending failed. turning it off.\n")
    APPEND = False

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

# Check MD suite
if(MD_PACKAGE == "amber"):
    from amber_module import MD_module
    md_module = MD_module(CONFIGFILE, DEBUG)
if(MD_PACKAGE == "gromacs"):
    print("Sorry, support for gromacs is not implemented yet.")
    sys.exit(-1)

# Initiate iterations
if APPEND:
    iterations = logger.loadIterations()
    #TODO: check if all files are present  
else:
    iterations.append(initiate.create_initial_iteration(INITIAL_TARGET_NUMBER_OF_SEGMENTS, 
                                                        INITIAL_BOUNDARIES, 
                                                        md_module, 
                                                        START_STATES, 
                                                        END_STATES))
    logger.log(iterations[0], CONFIGFILE)


#############################
#         Main Loop         #
#############################
for iteration_counter in range(len(iterations), MAX_ITERATIONS + 1):
    sys.stdout.write('\033[1m\nhdWE Status:\033[0m ' + 'Iteration '
    + str(iteration_counter).zfill(5) + '\n')    
    sys.stdout.flush()


    # 1. Sorting of segments into bins
    #    - initialize new iteration.
    #    - copy existing bins from the previous iteration into the new iteration.
    #    - Calculate the coordinate of each segment.
    #    - Generate new bins if required
    sys.stdout.write(' - Setting up new Iteration and sorting Segments into Bins\n')
    sys.stdout.flush()  
    
    iterations.append(Iteration(iteration_counter, iterations[-1].getBoundaries())) 
    resorting.copyBinStructureToLastIteration(iterations)
    resorting.resort(iterations, 
                     md_module,
                     INITIAL_TARGET_NUMBER_OF_SEGMENTS,
                     START_STATES,
                     END_STATES)
    
    # 2. Backup the segments lists of all bins
    #    - Saving the segment assignments for correct
    #      rate matrix calculation in the original_segments
    #      list of bins. 
    #    - This list should be immutable 
    for this_bin in iterations[-1]:
        this_bin.backupInitialSegments()
        
    # 3. Assign recycled probabilities
    if STEADY_STATE:
        recycling.recycleProbability(iterations[-1])
        print(" - Recyled Probability: {0:f}".format(iterations[-1].getProbabilityFlow()))
    
    # 4. Reweighting of bin probabilities
    #    The order of the following steps should no longer matter.  
    if iterations[-1].getNumberOfBins() > 1 and REWEIGHTING_RANGE > 0.0:
        sys.stdout.write(' - Reweighting Bin Probabilities\n')
        sys.stdout.flush()
        #TODO: why does reweighting need the workdir and jobname? fix it! 
        reweighting.reweightBinProbabilities(iterations,
                                             REWEIGHTING_RANGE,
                                             WORKDIR,
                                             JOBNAME)
        
    # 5. Resampling
    sys.stdout.write(' - Resampling\n')
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
    sys.stdout.write(' - Run MDs\n')
    sys.stdout.flush() 
    md_module.RunMDs(iterations[-1])
    
    
    # 6. log everything
    logger.log(iterations[-1], CONFIGFILE)

    #if DEBUG: 
    print(" - The overall probabiliy is {0:05f}".format(iterations[-1].getProbability()))

    #count total n of segments during iterations
    n_segments = 0
    for iteration_loop in iterations:
        n_segments += iteration_loop.getNumberOfPropagatedSegments()
    sys.stdout.write(' - (Total number of propagated segments: ' + str(n_segments-1)+ ')\n')
    sys.stdout.flush()
    
    #check for empty bins #TODO: make this a function of iterations
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
print('\nhdWE completed. Total number of propagated segments: ' + str(n_segments-1) + '         ')



