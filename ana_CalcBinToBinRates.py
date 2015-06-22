#!/usr/bin/env python2
"""
Return the bin to bin transition matrix or the the rate matrix
"""
from __future__ import print_function
import numpy
from lib.logger import Logger
import argparse
import lib.analysis_operations as analysis_operations

###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    metavar="FILE", required=True,
                    help="The log directory for reading")                  
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0,
                    help="First iteration to use.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Last iteration to use.")  
parser.add_argument('-t', '--bin-to-bin-transitions', action="store_true",
                    dest="bin_to_bin_transitions", 
                    help="Print a matrix with bin to bin transition events")
parser.add_argument('-o', '--output', dest="outfile", metavar="OUTPUT",
                    help="Write the matrix to this file")

args = parser.parse_args()
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration)

N = iterations[-1].getNumberOfBins()
transition_matrix = numpy.zeros((N,N), int)

if args.bin_to_bin_transitions:
    for iteration in iterations:
        for this_bin in iteration:
            for segment in this_bin.initial_segments:
                transition_matrix[iteration.bins[segment.getParentBinId()].getCoordinateIds(), 
                                  iteration.bins[segment.getBinId()].getCoordinateIds()]      += 1
    
    if not args.outfile: 
        print( transition_matrix )
    else:
        numpy.savetxt(args.outfile, transition_matrix, fmt='% 4d')
else:
    rateMatrix = analysis_operations.meanRateMatrix(iterations)
    if not args.outfile:
        print(rateMatrix)
    else:
        numpy.savetxt(args.outfile, rateMatrix, fmt='%.2f')
