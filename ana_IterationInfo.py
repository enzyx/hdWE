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

args = parser.parse_args()
logger = Logger(args.logdir)
try:
    iteration = logger.loadIteration(args.iteration_index)
except IOError:
    print("Could not find file for iteration {}!".format(args.iteration_index))
    sys.exit(-1)

for this_bin in iteration:
    print(" Bin: {:<4d} Segments: {:<4d} p: {:<5.4e}  Coord ids: {} Active: {}".format(this_bin.getId(),
                                               this_bin.getNumberOfSegments(),
                                               this_bin.getProbability(), 
                                               this_bin.getCoordinateIds(),
                                               this_bin.getSampleRegion()))