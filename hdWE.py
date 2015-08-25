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
import lib.reweighting as reweighting
import lib.cleanup as cleanup

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

# File Structure
WORKDIR                           = config.get('hdWE','WORKDIR')
JOBNAME                           = config.get('hdWE','JOBNAME')
LOGDIR                            = constants.getLogDirPath(WORKDIR, JOBNAME)
RUNDIR                            = constants.getRunDirPath(WORKDIR, JOBNAME)
# General Options
APPEND                            = args.append
APPEND_NEW_CONFIG                 = args.append_new_config
OVERWRITE                         = args.overwrite
DEBUG                             = args.debug
NUMBER_OF_THREADS                 = config_parser.parseNumberOfThreads(config)
KEEP_COORDS_FREQUENCY             = config_parser.parseKeepCoordsFrequency(config)
KEEP_COORDS_SEGMENTS              = config_parser.parseKeepCoordsSegments(config)
COMPRESS_ITERATION                = config_parser.parseCompressIteration(config)
COMPRESS_CLOSEST_MASK             = config_parser.parseCompressClosestMask(config)
# hdWE Options
STARTING_STRUCTURES               = config_parser.parseStartingStructures(config)
MAX_ITERATIONS                    = config_parser.parseMaxIterations(config)
INITIAL_TARGET_NUMBER_OF_SEGMENTS = int(config.get('hdWE','segments-per-bin'))
INITIAL_BOUNDARIES                = config_parser.parseInitialBoundaries(config)
INITIAL_SAMPLE_REGION             = config_parser.parseSampleRegion(config)
# -merging
MERGE_MODE                        = str(config.get('hdWE', 'merge-mode'))
if MERGE_MODE == 'closest':
    MERGE_THRESHOLD               = float(config.get('hdWE', 'merge-threshold'))
elif MERGE_MODE == 'weighted' or MERGE_MODE == 'random' or MERGE_MODE == 'marginonly':
    MERGE_THRESHOLD               = 0 # needs to be defined
else:
    raise Exception("Merge mode not found")    
# -reweighting
REWEIGHTING_RANGE                 = config_parser.parseReweightingRange(config)
if REWEIGHTING_RANGE > 0:
    REWEIGHTING_MAX_ITERATION     = config_parser.parseReweightingMaxIteration(config)

# MD module
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
    #sys.stderr.write("WARNING: appending failed. turning it off.\n")
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

# Initiate the reweighter
if REWEIGHTING_RANGE > 0:
    reweighter = reweighting.Reweighting( reweighting_range = REWEIGHTING_RANGE )

# Initiate iterations
if APPEND:
    iterations = logger.loadLastIterations(N=1)
    # Load the previous iterations to restore the rate matrices for the reweighter module
    if REWEIGHTING_RANGE > 0 and iterations[-1].getId() <= REWEIGHTING_MAX_ITERATION:
        print('Loading previous iterations to restore rate matrix for reweighting...')
        for iteration_counter_tmp in range(1,iterations[-1].getId() + 1):
            iteration_tmp = logger.loadIteration(iteration_counter_tmp)
            reweighter.storeRateMatrix(iteration_tmp)
        iteration_tmp = []
            
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

# Handle the deletion/compression of MD output files 
cleaner = cleanup.Cleanup(md_module, NUMBER_OF_THREADS, COMPRESS_ITERATION, 
                 COMPRESS_CLOSEST_MASK, KEEP_COORDS_FREQUENCY, KEEP_COORDS_SEGMENTS, DEBUG)

#############################
#         Main Loop         #
#############################
for iteration_counter in range(iterations[-1].getId() + 1, MAX_ITERATIONS + 1):
    sys.stdout.write('\033[1m\nhdWE\033[0m Job: {jn} Status: Iteration {it:05d}\n'.format(jn = JOBNAME, it = iteration_counter))
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
    
    if MERGE_MODE == "closest":
        bins_rmsds = md_module.calcBinsRmsds(iterations[-1])
        for this_bin in iterations[-1]:
            this_bin.resampleSegments(MERGE_MODE, 
                                      MERGE_THRESHOLD, 
                                      len(iterations[-1].boundaries[0]),
                                      bins_rmsds[this_bin.getId()])
    else:
        # Parallel
        if NUMBER_OF_THREADS > 1:
            thread_container = ThreadContainer()
            for this_bin in iterations[-1]:
                thread_container.appendJob(threading.Thread(target=this_bin.resampleSegments, args=(MERGE_MODE, 
                                                                                             MERGE_THRESHOLD,
                                                                                             len(iterations[-1].boundaries[0]) )))
                if thread_container.getNumberOfJobs() >= NUMBER_OF_THREADS:
                    thread_container.runJobs()
            # Run remaining jobs
            thread_container.runJobs()
        # Serial
        else:
            for this_bin in iterations[-1]:
                this_bin.resampleSegments(MERGE_MODE, 
                                          MERGE_THRESHOLD,
                                          len(iterations[-1].boundaries[0]))
  
    # 4. Treatment of outer-region bins.
    #    Has to be after resampling for consistency with ana_TraceFlux  
    iterations[-1].resetOuterRegion()          

    # 5. Reweighting
    if REWEIGHTING_RANGE > 0.0 and iteration_counter <= REWEIGHTING_MAX_ITERATION:
        sys.stdout.write(' - Reweighting\n')
        reweighter.storeRateMatrix(iterations[-1])
        if iterations[-1].getNumberOfBins() > 1:
            reweighter.reweightBinProbabilities(iterations[-1])
                
    # 6. Run MDs
    sys.stdout.write(' - Run MDs\n')
    sys.stdout.flush() 
    md_module.runMDs(iterations[-1])
    sys.stdout.write('\n')
    sys.stdout.flush() 
    
    # 7. Calculate Segment Coordinates
    sys.stdout.write(' - Calculate Coordinates\n')
    sys.stdout.flush() 
    md_module.calcCoordinates(iterations[-1])    
    
    # 8. log everything
    logger.log(iterations[-1], CONFIGFILE)

    # 9. compress files and delete unwanted files
    print(" - Compressing/Deleting MD files")    
    if len(iterations) > 2:
        cleaner.doCleanup(iterations[-3])
    
    if DEBUG: 
        print("\n    The overall probability is {0:05f}".format(iterations[-1].getProbability()))
    
    print('    Empty bins: ' + str(iterations[-1].countEmptyBins()))
    
    # 10. Save some RAM
    iterations = iterations[-3:]

#############################
#     End of Main Loop      #
#############################
cleaner.finishCleanupThreads()
print('\nhdWE completed.')



