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
prev_data = parser.add_mutually_exclusive_group()
prev_data.add_argument('-s', '--segments', dest="print_segments", default=False, action="store_true",
                    help="Print segment info")
prev_data.add_argument('-p', '--parent-segments', dest="parent_segments", default=False, action="store_true",
                    help="Print segment parent info")


args = parser.parse_args()
logger = Logger(args.logdir)
try:
    iteration = logger.loadIteration(args.iteration_index)
except IOError:
    print("Could not find file for iteration {}!".format(args.iteration_index))
    sys.exit(-1)

def bin_compare(x, y):
    for dimension in range(0, len(x.coordinate_ids)):
        diff = (x.coordinate_ids[dimension] - y.coordinate_ids[dimension])
        if (x.coordinate_ids[dimension] - y.coordinate_ids[dimension]) != 0:
            return diff

#iteration.bins = sorted(iteration.bins, cmp = bin_compare)


for this_bin in iteration:
    print(" Bin: {:<4d} Segments: {:<4d} p: {:<5.4e}  Coord ids: {} Active: {}".format(this_bin.getId(),
                                               this_bin.getNumberOfSegments(),
                                               this_bin.getProbability(), 
                                               this_bin.getCoordinateIds(),
                                               this_bin.getSampleRegion()))
    if args.print_segments and not args.parent_segments:
        for segment in this_bin:
            print("    Segment: {:<4d} p: {:<5.4e} Coords: {}".format(segment.getId(),
                                                                      segment.getProbability(),
                                                                      segment.getCoordinates()))
    elif args.parent_segments:
        for segment in this_bin:
            print("    Segment: {:<4d} Parent bin id: {:<4d} Parent segment id: {:<4d}".format(segment.getId(),
                                                                                 segment.getParentBinId(),
                                                                                 segment.getParentSegmentId()))