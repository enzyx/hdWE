#!/usr/bin/python2
"""
hdWE is a hyperdimensional weighted ensemble simulation implementation.
"""
from __future__ import print_function
import sys
import ConfigParser
import argparse
import threading
import os
import glob
from   lib.thread_container import ThreadContainer
from   lib.iteration import Iteration
from   lib.logger import Logger 
import lib.initiate as initiate
import lib.constants as constants
import lib.resorting as resorting
import lib.recycling as recycling

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
NUMBER_OF_THREADS                 = int(config.get('hdWE','number-of-threads'))
KEEP_COORDS_FREQUENCY             = int(config.get('hdWE', 'keep-coords-frequency'))

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
    from lib.amber_module import MD_module
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
                                                        md_module))
    logger.log(iterations[0], CONFIGFILE)


#############################
#         Main Loop         #
#############################
for iteration_counter in range(iterations[-1].getId() + 1, MAX_ITERATIONS + 1):
    sys.stdout.write('\033[1m\nhdWE Status:\033[0m Iteration {:05d}\n'.format(iteration_counter))
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
                     INITIAL_TARGET_NUMBER_OF_SEGMENTS)
    
    # 2. Backup the segments lists of all bins
    #    - Saving the segment assignments for correct
    #      rate matrix calculation in the original_segments
    #      list of bins. 
    #    - This list should be immutable 
    for this_bin in iterations[-1]:
        this_bin.backupInitialSegments()

    # 3. Resampling
    sys.stdout.write(' - Resampling\n')
    sys.stdout.flush() 
    # Parallel
    if NUMBER_OF_THREADS > 1:
        thread_container = ThreadContainer()
        for this_bin in iterations[-1]:
            thread_container.appendJob(threading.Thread(target=this_bin.resampleSegments() ))
            if thread_container.getNumberOfJobs() >= NUMBER_OF_THREADS:
                thread_container.runJobs()
        # Run remaining jobs
        thread_container.runJobs()
    # Serial
    else:
        for this_bin in iterations[-1]:
            this_bin.resampleSegments()

    # 4. Run MDs
    sys.stdout.write(' - Run MDs\n')
    sys.stdout.flush() 
    md_module.RunMDs(iterations[-1])
    sys.stdout.write('\n')
    sys.stdout.flush() 
    
    # 5. log everything
    logger.log(iterations[-1], CONFIGFILE)

    # 6. delete unwanted files
    print(" - Deleting md files")
    if iterations[-2].getId() % KEEP_COORDS_FREQUENCY != 0:
        md_module.removeCoordinateFiles(iterations[-2])

    #if DEBUG: 
    print("\n    The overall probability is {0:05f}".format(iterations[-1].getProbability()))
    
    #check for empty bins #TODO: make this a function of iterations
    empty_bins = 0
    for bin_loop in iterations[-1]:
        if bin_loop.getNumberOfSegments() == 0:
            empty_bins += 1
    print('    Empty bins: ' + str(empty_bins))
    
    # Save some RAM
    iterations = iterations[-2:]

#############################
#     End of Main Loop      #
#############################
print('\nhdWE completed.')



