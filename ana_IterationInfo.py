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
                    help="Iteration which is displayed")
parser.add_argument('-n', '--initial-segments', dest="initial_segments", default=False, action="store_true",
                    help="Print initial segments info")
prev_data = parser.add_mutually_exclusive_group()
prev_data.add_argument('-s', '--segments', dest="print_segments", default=False, action="store_true",
                    help="Print segment info")
prev_data.add_argument('-p', '--parent-segments', dest="parent_segments", default=False, action="store_true",
                    help="Print segment parent info")
prev_data.add_argument('--show-inactive-bins', dest="show_inactive_bins", default=False, action="store_true",
                    help="Show data also for inactive bins")

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

iteration.bins = sorted(iteration.bins, cmp = bin_compare)

for this_bin in iteration:
    probs = []
    for this_segment in this_bin:
        probs.append(this_segment.getProbability())
    if len(probs)> 0:
        this_bin.probability_range = 1.0 *max(probs) / min(probs)
    else:
        this_bin.probability_range = 0


#######################################
#      Print initial segments         #
#######################################
if args.initial_segments:
    for this_bin in iteration:
        if args.show_inactive_bins == False and this_bin.getSampleRegion() == False:
            continue
        print(" Bin: {:<4d} Segments: {:<4d} p: {:<5.4e}  Coord ids: {} Active: {} Prob range: {:<3.2e}".format(this_bin.getId(),
                                                   this_bin.getNumberOfInitialSegments(),
                                                   this_bin.getInitialProbability(), 
                                                   this_bin.getCoordinateIds(),
                                                   this_bin.getSampleRegion(),
                                                   this_bin.probability_range))
        if args.print_segments:
            for segment in this_bin.initial_segments:
                    print("    Segment: {:<4d} p: {:<5.4e} Coords: {}".format(segment.getId(),
                                                                              segment.getProbability(),
                                                                              segment.getCoordinates()))
    sys.exit()
#######################################
for this_bin in iteration:
    if args.show_inactive_bins == False and this_bin.getSampleRegion() == False:
        continue
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