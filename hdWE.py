#!/usr/bin/python2
"""
  _         ___          ________ 
 | |       | \ \        / /  ____|
 | |__   __| |\ \  /\  / /| |__   
 | '_ \ / _` | \ \/  \/ / |  __|  
 | | | | (_| |  \  /\  /  | |____ 
 |_| |_|\__,_|   \/  \/   |______|
                                  
A hyperdimensional weighted ensemble simulation implementation.
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
import lib.config_parser as config_parser
import lib.constants as constants
import lib.resorting as resorting

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
                       help='continue previous iterations (use parameters from conf file when --read is given)')
prev_data.add_argument('-o', '--overwrite', dest="overwrite", action='store_true', 
                       default=False,
                       help='overwrite previous simulation data')
parser.add_argument('-n','--new-conf', dest="append_new_config", action='store_true', 
                       default=False,
                       help='Read new boundaries from config file when --append is used')

args = parser.parse_args()
if args.append_new_config and not args.append:
    parser.error('--new-conf can only be set together with --append.')

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
APPEND_NEW_CONFIG                 = args.append_new_config
OVERWRITE                         = args.overwrite
DEBUG                             = args.debug
MAX_ITERATIONS                    = int(config.get('hdWE','max-iterations'))
INITIAL_TARGET_NUMBER_OF_SEGMENTS = int(config.get('hdWE','segments-per-bin'))
INITIAL_BOUNDARIES                = config_parser.parseInitialBoundaries(config)
INITIAL_SAMPLE_REGION             = config_parser.parseSampleRegion(config)
STARTING_STRUCTURES               = config_parser.parseStartingStructures(config)
NUMBER_OF_THREADS                 = int(config.get('hdWE','number-of-threads'))
KEEP_COORDS_FREQUENCY             = config_parser.parseKeepCoordsFrequency(config)
KEEP_COORDS_SEGMENTS              = config_parser.parseKeepCoordsSegments(config)
MERGE_MODE                        = str(config.get('hdWE', 'merge-mode'))
MERGE_THRESHOLD                   = float(config.get('hdWE', 'merge-threshold'))

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
    iterations = logger.loadLastIterations(N=1)
    if APPEND_NEW_CONFIG:
        iterations[-1].boundaries    = INITIAL_BOUNDARIES
        iterations[-1].sample_region = INITIAL_SAMPLE_REGION
        iterations[-1].target_number_of_segments  = INITIAL_TARGET_NUMBER_OF_SEGMENTS
        for this_bin in iterations[-1]:
            this_bin.target_number_of_segments = iterations[-1].target_number_of_segments
            this_bin.sample_region = iterations[-1].isInSampleRegion(this_bin.getCoordinateIds()) 

    #TODO: check if all files are present
else:
    iterations.append(initiate.createInitialIteration(STARTING_STRUCTURES,
                                                      WORKDIR, 
                                                      JOBNAME,
                                                      INITIAL_TARGET_NUMBER_OF_SEGMENTS, 
                                                      INITIAL_BOUNDARIES, 
                                                      INITIAL_SAMPLE_REGION,
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
    
    iterations.append(Iteration(iteration_counter, 
                                iterations[-1].getBoundaries(), 
                                iterations[-1].getSampleRegion()) ) 
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
            thread_container.appendJob(threading.Thread(target=this_bin.resampleSegments(MERGE_MODE, 
                                                                                         MERGE_THRESHOLD, 
                                                                                         md_module) ))
            if thread_container.getNumberOfJobs() >= NUMBER_OF_THREADS:
                thread_container.runJobs()
        # Run remaining jobs
        thread_container.runJobs()
    # Serial
    else:
        for this_bin in iterations[-1]:
            this_bin.resampleSegments(MERGE_MODE, 
                                      MERGE_THRESHOLD, 
                                      md_module)
  
    # 4. Treatment of outer-region bins.
    #    Has to be after resampling for consistency with ana_TraceFlux  
    iterations[-1].resetOuterRegion()          
            

    # 5. Run MDs
    sys.stdout.write(' - Run MDs\n')
    sys.stdout.flush() 
    md_module.runMDs(iterations[-1])
    sys.stdout.write('\n')
    sys.stdout.flush() 
    
    # 6. Calculate Segment Coordinates
    sys.stdout.write(' - Calculate Coordinates\n')
    sys.stdout.flush() 
    md_module.calcCoordinates(iterations[-1])    
    
    # 7. log everything
    logger.log(iterations[-1], CONFIGFILE)

    # 8. delete unwanted files
    print(" - Deleting md files")
    if iterations[-2].getId() % KEEP_COORDS_FREQUENCY != 0:
        md_module.removeCoordinateFiles(iterations[-2])
    else:
        if KEEP_COORDS_SEGMENTS > 0:
            md_module.removeCoordinateFiles(iterations[-2], KEEP_COORDS_SEGMENTS)
    
    if DEBUG: 
        print("\n    The overall probability is {0:05f}".format(iterations[-1].getProbability()))
    
    #check for empty bins #TODO: make this a function of iterations
    empty_bins = 0
    for bin_loop in iterations[-1]:
        if bin_loop.getNumberOfSegments() == 0 and bin_loop.getSampleRegion == True:
            empty_bins += 1
    print('    Empty bins: ' + str(empty_bins))
    
    # Save some RAM
    iterations = iterations[-2:]

#############################
#     End of Main Loop      #
#############################
print('\nhdWE completed.')



