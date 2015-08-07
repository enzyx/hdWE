#!/usr/bin/env python2
"""
Return the bin to bin transition matrix or the the rate matrix
"""
from __future__ import print_function
import numpy
from lib.logger import Logger
import argparse
import lib.analysis_operations as analysis_operations
import sys

###### Function ######
def printMatrix(matrix):
    digits = str(len(str(matrix.max())) + 1)
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if i==j:
                sformat = "\033[0;31m{: "+digits+"d}\033[0m"
            elif i==j+1 or i==j-1:
                sformat = "\033[0;32m{: "+digits+"d}\033[0m"
            else:
                sformat = "{: "+digits+"d}"
            sys.stdout.write(sformat.format(matrix[i][j]))
        sys.stdout.write('\n')

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
parser.add_argument('-s', '--sort', dest='sort_dimension', 
                    metavar='INT', default=0, type=int, 
                    help='Dimension index along which to sort bins.')

args = parser.parse_args()
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration)

N = iterations[-1].getNumberOfBins()
DIMENSION = args.sort_dimension
transition_matrix = numpy.zeros((N,N), int)

if args.bin_to_bin_transitions:
    for iteration in iterations:
        for this_bin in iteration:
            for segment in this_bin.initial_segments:
                print(segment.getParentBinId(), iteration.getNumberOfBins())
                transition_matrix[iteration.bins[segment.getParentBinId()].getId(), 
                                  iteration.bins[segment.getBinId()].getId()]      += 1

    if args.sort_dimension != -1:                       
        # sort along dimension
        sort_matrix = []
        largest_coordinate_id = 0
        for this_bin in iterations[-1]:
            if this_bin.getCoordinateIds()[DIMENSION] > largest_coordinate_id:
                largest_coordinate_id = this_bin.getCoordinateIds()[DIMENSION]
        sort_dimensions = numpy.max([(largest_coordinate_id+1), N])
        large_sort_matrix = numpy.zeros((sort_dimensions, N), int)
        for this_bin in iterations[-1]:
            this_coord = this_bin.getCoordinateIds()[DIMENSION]
            large_sort_matrix[this_coord, this_bin.getId()] = 1
        b_got_values = False
        # cut all empty bins before first sampled one, 
        for line in large_sort_matrix:
            if numpy.sum(line) != 0:
                b_got_values = True
            if not b_got_values:
                continue    
            sort_matrix.append(line)
        sort_matrix = numpy.asarray(sort_matrix, dtype=int)
        # find lines full of zeros and remember their index
        unfound_coord_ids = []
        for line_index, line in enumerate(sort_matrix):
            if numpy.sum(line) == 0:
                unfound_coord_ids.append(line_index)
    
         
        transition_matrix = numpy.dot(sort_matrix, transition_matrix)
        transition_matrix = numpy.dot(transition_matrix, numpy.transpose(sort_matrix))
    
        # set values of unfound bins to -1
        for coord_id in unfound_coord_ids:
            transition_matrix[coord_id, :] = -1

    if not args.outfile: 
        printMatrix( transition_matrix )
    else:
        numpy.savetxt(args.outfile, transition_matrix, fmt='% 4d')
else:
    rateMatrix = analysis_operations.meanRateMatrix(iterations)
    if not args.outfile:
        print(rateMatrix)
    else:
        numpy.savetxt(args.outfile, rateMatrix, fmt='%.2f')

# start_start_rate = 0.0
# for i in range(len(rateMatrix)):
#     for j in range(len(rateMatrix)):
#         if iterations[-1].bins[j].isStartStateBin() == True:
#             start_start_rate += rateMatrix[i,j] * iterations[-1].bins[i].getProbability()
# 
# print (start_start_rate)
        

