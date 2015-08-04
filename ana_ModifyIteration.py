#!/usr/bin/env python2
"""
Modify an existing simulation by skipping initial bins
"""
from __future__ import print_function
from lib.logger import Logger
import argparse
import sys

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    metavar="FOLDER", required=True,
                    help="The log directory for reading")
parser.add_argument('-i', '--iteration', dest="iteration_index",
                    type=int, default=-1,
                    help="Iteration which is modified")  
parser.add_argument('-t', '--target-number-of-bins', dest="bins_target", nargs='+', 
                    metavar="BIN TARGET", type=int,
                    required=False, help="Change target number of segments for bins in format [binId, targetNumerOfSegments, ...]")
parser.add_argument('-s', '--skip', dest="skip_bins",
                    type=int, default=0, 
                    help="Number of bins to be skipped upon copying.")
parser.add_argument('-f', '--overwrite', dest="overwrite",
                    default=False,  action='store_true',
                    help="Overwrite output files.")


args = parser.parse_args()
logger = Logger(args.logdir)
try:
    iteration = logger.loadIteration(args.iteration_index)
except IOError:
    print("Could not find file for iteration {}!".format(args.iteration_index))
    sys.exit(-1)

if not args.overwrite:
    print("Running dry! Use -f to save changes.")

# Change target number of segments if required
if args.bins_target:
    for binId, targetNumberOfSegments in zip(*[iter(args.bins_target)]*2):
        print("bin: {} targetNumberOfSegments: {} --> {}".format(binId, 
                                  iteration.bins[binId].getTargetNumberOfSegments(), 
                                  targetNumberOfSegments))
        iteration.bins[binId].target_number_of_segments = targetNumberOfSegments

if args.overwrite:
    print("Overwriting iteration file: ", args.iteration_index)
    logger.logIteration(iteration)