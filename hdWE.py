#!/usr/bin/python2
#
# This file is part of hdWE. 
# Copyright (C) 2016 Manuel Luitz <manuel.luitz@tum.de>
# Copyright (C) 2016 Rainer Bomblies <r.bomblies@tum.de>
# Copyright (C) 2016 Fabian Zeller
#
# hdWE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# hdWE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with hdWE. If not, see <http://www.gnu.org/licenses/>.
# 
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
import os
import glob
from   lib.iteration import Iteration
from   lib.logger import Logger 
import lib.initiate as initiate
import lib.config_parser as config_parser
import lib.constants as constants
import lib.resorting as resorting
import lib.reweighting as reweighting
import lib.cleanup as cleanup
import lib.resampling as resampling

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
CALCULATE_VELOCITIES              = config_parser.parseCalculateVelocities(config) 
# hdWE Options
STARTING_STRUCTURES               = config_parser.parseStartingStructures(config)
MAX_ITERATIONS                    = config_parser.parseMaxIterations(config)
INITIAL_TARGET_NUMBER_OF_SEGMENTS = config_parser.parserInitialNumberOfTargetSegments(config)
INITIAL_BOUNDARIES                = config_parser.parseInitialBoundaries(config)
INITIAL_SAMPLE_REGION             = config_parser.parseSampleRegion(config)
STEADY_STATE                      = config_parser.parseSteadyState(config)
START_BIN_COORDINATE_IDS          = config_parser.parseStartBinCoordinateIds(config)
# resampling
RESAMPLING_MODE                   = config_parser.parseResamplingMode(config)
CLOSEST_MERGE_THRESHOLD           = config_parser.parseMergeThreshold(config)
SPLIT_FORWARD_NUMBER_OF_CHILDREN  = config_parser.parseSplitForwardNumberOfChildren(config)
PRIMARY_COORDINATE                = config_parser.parsePrimaryCoordinate(config)
SPLIT_REGION                      = config_parser.parseSplitRegion(config)
FRONT_INTERVAL                    = config_parser.parseFrontInterval(config)
# reweighting
REWEIGHTING_RANGE                 = config_parser.parseReweightingRange(config)
REWEIGHTING_MAX_ITERATION         = config_parser.parseReweightingMaxIteration(config)

# MD module
if "amber" in config.sections():
    MD_PACKAGE = "amber"
elif "gromacs" in config.sections():
    MD_PACKAGE = "gromacs"
elif "langevin" in config.sections():
    MD_PACKAGE = "langevin"
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
if(MD_PACKAGE == "langevin"):
    from lib.langevin_module import MD_module
    md_module = MD_module(CONFIGFILE, DEBUG)
if(MD_PACKAGE == "gromacs"):
    print("Sorry, support for gromacs is not implemented yet.")
    sys.exit(-1)

# Initiate the reweighter
if REWEIGHTING_RANGE > 0:
    reweighter = reweighting.Reweighting( reweighting_range = REWEIGHTING_RANGE )

# Initiate iterations
if APPEND == True:
    # load the last two iterations if existing for proper function of
    # cleanup module
    if logger.getLastIterationId() > 1: 
        iterations = logger.loadLastIterations(N=2)
    else:
        logger.loadLastIterations(N=1)
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
        for this_bin in iterations[-1].bins:
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
                                                      md_module,
                                                      START_BIN_COORDINATE_IDS))
    logger.log(iterations[0], CONFIGFILE)  

# Create an instance of the resampling module
resampler = resampling.Resampling(md_module, RESAMPLING_MODE, CLOSEST_MERGE_THRESHOLD, 
                                  PRIMARY_COORDINATE, SPLIT_FORWARD_NUMBER_OF_CHILDREN, 
                                  SPLIT_REGION, FRONT_INTERVAL)

# Handle the deletion/compression of MD output files 
cleaner = cleanup.Cleanup(md_module, NUMBER_OF_THREADS, COMPRESS_ITERATION, 
                          COMPRESS_CLOSEST_MASK, KEEP_COORDS_FREQUENCY, 
                          KEEP_COORDS_SEGMENTS, DEBUG)

#############################
#         Main Loop         #
#############################
for iteration_counter in range(iterations[-1].getId() + 1, MAX_ITERATIONS + 1):
    sys.stdout.write('\033[1m\nhdWE\033[0m Job: {jn} Status: Iteration {it:08d}\n'.format(jn = JOBNAME, it = iteration_counter))
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
                                iterations[-1].getSampleRegion(),
                                iterations[-1].getNumberOfStartingStructures() )) 
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
    resampler.resample(iterations[-1])                                     
  
    # 4. Treatment of outer-region bins.
    #    Has to be after resampling for consistency with ana_TraceFlux  
    iterations[-1].resetOuterRegion(STEADY_STATE)          

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
    if CALCULATE_VELOCITIES:
        sys.stdout.write(' - Calculate Velocities\n')
        md_module.calcVelocity(iterations[-1])   
    
    # 8. log everything
    logger.log(iterations[-1], CONFIGFILE)

    # 9. compress files and delete unwanted files
    print(" - Cleanup of MD files")    
    if len(iterations) > 2:
        cleaner.doCleanup(iterations[-3])
    
    if DEBUG: 
        print("\n    The overall probability is {0:05f}".format(iterations[-1].getProbability()))
    
    # 10. Save some RAM
    iterations = iterations[-3:]

#############################
#     End of Main Loop      #
#############################
cleaner.finishCleanupThreads()
print('\nhdWE completed.')



