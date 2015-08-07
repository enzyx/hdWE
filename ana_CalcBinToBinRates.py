#!/usr/bin/env python2
"""
Return the bin to bin transition matrix or the the rate matrix
"""
from __future__ import print_function
import numpy as np
from lib.logger import Logger
import argparse
import lib.analysis_operations as analysis_operations
import sys

###### Function ######

def sortMatrix(matrix, order):
    """
    sorts rows and columns of a matrix according
    to given order. order is interpreted as order[source_index] = target_index
    """
    # if order has wrong dimension or 
    if len(order) != len(matrix):
        raise Exception('wrong order dimension given')
    # or matrix isn't square I can't work
    elif not all (len (row) == len (matrix) for row in matrix):
        raise Exception('can\'t sort non-square matrix')
    
    sort_matrix = np.zeros((len(order),len(order)), dtype=int)
    for source_index, target_index in enumerate(order):
        sort_matrix[target_index, source_index] = 1

         
    new_matrix = np.dot(sort_matrix, matrix)
    new_matrix = np.dot(new_matrix, np.transpose(sort_matrix))        
    return new_matrix
        

def printMatrix(matrix, coord_ids = []):
    digits = str(len(str(matrix.max())) + 1)
    if len(coord_ids) != 0:
        sys.stdout.write("coordinate id of |   transition\nsort dimension   |     matrix\n")
    for i in range(len(matrix)):
        if len(coord_ids) != 0:
            sys.stdout.write("{:2d} ".format(int(coord_ids[i])))
            sys.stdout.write("|")
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
                    help='Dimension index along which to sort bins. -1 does not sort.')

args = parser.parse_args()
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration)

N = iterations[-1].getNumberOfBins()
DIMENSIONS = len(iterations[-1].bins[0].getCoordinateIds())
SORT_DIMENSION = args.sort_dimension
transition_matrix = np.zeros((N,N), int)

if args.bin_to_bin_transitions:
    for iteration in iterations:
        for this_bin in iteration:
            for segment in this_bin.initial_segments:
                transition_matrix[iteration.bins[segment.getParentBinId()].getId(), 
                                  iteration.bins[segment.getBinId()].getId()]      += 1

    if args.sort_dimension != -1:    
        # collect coordinate ids of bins        
        bin_list = []
        for this_bin in iteration:
            # save bin_id as additional element to know sorting order later
            coord_id = this_bin.getCoordinateIds()[:].tolist()
            coord_id.append(this_bin.getId())
            bin_list.append(np.asarray(coord_id))

        # sort in order of dimensions
        sort_dimension_order = [SORT_DIMENSION] \
                               + np.arange(SORT_DIMENSION).tolist() \
                               + np.arange(SORT_DIMENSION+1, DIMENSIONS).tolist()
        sorted_bin_list = np.asarray(sorted(bin_list, 
                                            key=lambda element: tuple([element[i] for i in sort_dimension_order])))
        sorting_order = sorted_bin_list[:,-1]
        
        transition_matrix = sortMatrix(transition_matrix, sorting_order)
        coordinate_ids_for_plot = sorted_bin_list[:,SORT_DIMENSION]
    else:
        # to preserve non-sorting functionality          
        coordinate_ids_for_plot = []
    if not args.outfile: 
        printMatrix( transition_matrix, coordinate_ids_for_plot)
    else:
        np.savetxt(args.outfile, transition_matrix, fmt='% 4d')
else:
    rateMatrix = analysis_operations.meanRateMatrix(iterations)
    if not args.outfile:
        print(rateMatrix)
    else:
        np.savetxt(args.outfile, rateMatrix, fmt='%.2f')


