#!/usr/bin/env python2
"""
Modify an existing simulation by skipping initial bins
"""
from __future__ import print_function
import sys
import os
import shutil
from logger import Logger
import argparse  
from  iteration import Iteration
import constants
import ConfigParser

###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    metavar="FOLDER", required=True,
                    help="The log directory for reading")
parser.add_argument('-r', '--run', type=str, dest="rundir", 
                    metavar="FOLDER", required=True,
                    help="The run directory for reading")      
parser.add_argument('-i', '--iteration', dest="iteration_index",
                    type=int, default=-1,
                    help="Iteration which is copied and modified")  
parser.add_argument('-j', '--jobname', dest="mod_jobname", 
                    type=str, required=True,
                    help="Output jobname")
parser.add_argument('-s', '--skip', dest="skip_bins",
                    type=int, default=0, 
                    help="Number of bins to be skipped upon copying.")
parser.add_argument('-f', '--overwrite', dest="overwrite",
                    default=False,  action='store_true',
                    help="Overwrite output files.")
parser.add_argument('-u', '--uniform', dest="uniform",
                    default=False,  action='store_true',
                    help="Uniform probability distribution")

args = parser.parse_args()
WORKDIR         = os.getcwd() + "/"
NEW_JOBNAME     = args.mod_jobname
NEW_LOGDIR      = constants.getLogDirPath(WORKDIR, NEW_JOBNAME)
NEW_RUNDIR      = constants.getRunDirPath(WORKDIR, NEW_JOBNAME)
OLD_LOGDIR      = args.logdir + "/"
OLD_RUNDIR      = args.rundir + "/"
ITERATION_INDEX = args.iteration_index
SKIP_BINS       = args.skip_bins
UNIFORM         = args.uniform

#get the actual Iteration from logger module
logger = Logger(OLD_LOGDIR)
old_iteration    = logger.loadIteration(ITERATION_INDEX)
OLD_CONFIG_FNAME = logger.loadConfigFile(ITERATION_INDEX)
NEW_CONFIG_FNAME = WORKDIR + "/" + NEW_JOBNAME + ".conf"

# Load the config file
old_config = ConfigParser.ConfigParser()
old_config.read(OLD_CONFIG_FNAME)

# Print some interesting informations
print("Found config file: {}".format(OLD_CONFIG_FNAME))
print("Using iteration {0}".format(old_iteration.getId()))
print("Found {0} bins.".format(old_iteration.getNumberOfBins()))
print("Saving everything to new jobname: ", NEW_JOBNAME)

# Check files and folders
if not args.overwrite:
    if os.path.exists(NEW_JOBNAME + ".conf"):
        print("Configuration file already exists: {0}".format(NEW_JOBNAME + ".conf" ))
        print("Use --overwrite")
        sys.exit()
    if os.path.exists(NEW_LOGDIR):
        print("Log folder exists: {0}".format(NEW_LOGDIR))
        print("Use --overwrite")
        sys.exit()
    if os.path.exists(NEW_RUNDIR):
        print("Run folder exists: {0}".format(NEW_RUNDIR))
        print("Use --overwrite")
        sys.exit()

try: os.remove(NEW_JOBNAME + ".conf") 
except OSError: pass
try: os.remove(NEW_JOBNAME + ".log") 
except OSError: pass
shutil.rmtree(NEW_RUNDIR, ignore_errors=True) 
shutil.rmtree(NEW_LOGDIR, ignore_errors=True) 

os.mkdir(NEW_RUNDIR)
os.mkdir(NEW_LOGDIR)

# Generate a new iteration
iteration = Iteration(0)
bin_index = 0
for _bin in old_iteration:
    if _bin.getId() < SKIP_BINS:
        continue
    iteration.generateBin(reference_iteration_id=0, 
                          reference_bin_id=bin_index, 
                          reference_segment_id=0, 
                          target_number_of_segments=_bin.getTargetNumberOfSegments(), 
                          outrates_converged=False)
    iteration.bins[-1].generateSegment(probability=_bin.getProbability(),
                                       parent_iteration_id = 0,
                                       parent_bin_id=bin_index,
                                       parent_segment_id=0)
    
    print(_bin.getReferenceNameString())
    shutil.copy2(OLD_RUNDIR + _bin.getReferenceNameString() + ".rst7",
                 NEW_RUNDIR + iteration.bins[-1].getReferenceNameString() + ".rst7")
    _bin.backupInitialSegments()
    bin_index += 1

# Fix probabilities
p_tot = iteration.getProbability()
for _bin in iteration:
    for segment in _bin:
        segment.setProbability(segment.getProbability() / p_tot)

if UNIFORM:
    p_bin = 1.0/iteration.getNumberOfBins()
    for _bin in iteration:
        for segment in _bin:
            segment.setProbability(p_bin)

# Save everything
new_logger = Logger(NEW_LOGDIR)
new_logger.logIteration(iteration)

# Write hdWE config file
new_config = ConfigParser.RawConfigParser()
new_config.read(OLD_CONFIG_FNAME)
new_config.set('hdWE', 'workdir' , WORKDIR)
new_config.set('hdWE', 'jobname', NEW_JOBNAME)
with open(NEW_CONFIG_FNAME, 'wb') as CFILE:
    new_config.write(CFILE)
CFILE.close() 