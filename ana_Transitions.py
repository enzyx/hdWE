#!/usr/bin/env python2
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
Return the bin to bin transition matrix or the the rate matrix
"""
from __future__ import print_function
import numpy as np
from lib.logger import Logger
import argparse
import lib.analysis_operations as analysis_operations
import sys

###### Function ######

def sortList(this_list, order):
    """
    sorts a list according to given order
    order is interpreted as order[source_index] = target_index
    """
    new_list = [0] * len(this_list)
    for old_index, new_index in enumerate(order):
        new_list[new_index] = list[old_index] 
    return new_list

def sortMatrix(matrix, order):
    """
    sorts rows and columns of a matrix according
    to given order. order is interpreted as order[target_index] = source_index
    """
    # if order has wrong dimension or 
    if len(order) != len(matrix):
        raise Exception('wrong order dimension given')
    # or matrix isn't square I can't work
    elif not all (len (row) == len (matrix) for row in matrix):
        raise Exception('can\'t sort non-square matrix')
    
    sort_matrix = np.zeros((len(order),len(order)), dtype=int)
    for target_index, source_index in enumerate(order):
        sort_matrix[target_index, source_index] = 1
                        
    new_matrix = np.dot(sort_matrix, matrix)
    new_matrix = np.dot(new_matrix, np.transpose(sort_matrix))        
    return new_matrix
        

def printMatrix(matrix, coord_ids = []):
    digits = str(len(str(matrix.max())) + 1)
    # coord id header
    if len(coord_ids) != 0:
        sys.stdout.write("coordinate ids | transition matrix\n")
    for i in range(len(matrix)):
        # coordinate ids:
        if len(coord_ids) != 0:
            sys.stdout.write("{:14s} ".format(str(coord_ids[i])))
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
parser.add_argument('-s', '--sort', nargs='?', type=int,
                    metavar='INT', default=False, const=0, 
                    help='Dimension index along which to sort bins.')

args = parser.parse_args()
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration, verbose=True)

N_BINS = iterations[-1].getNumberOfBins()
DIMENSIONS = len(iterations[-1].bins[0].getCoordinateIds())
transition_matrix = np.zeros((N_BINS,N_BINS), int)

if args.bin_to_bin_transitions:
    for iteration in iterations:
        for this_bin in iteration:
            for segment in this_bin.initial_segments:
                transition_matrix[iteration.bins[segment.getParentBinId()].getId(), 
                                  iteration.bins[segment.getBinId()].getId()]      += 1

    # to preserve non-sorting functionality  
    coordinate_ids_for_plot = []
    
    for this_bin in iterations[-1]:
        coordinate_ids_for_plot.append(this_bin.getCoordinateIds())

    # if must be so clumsy because python evaluates 0 as false
    if type(args.sort) == int:
        SORT_DIMENSION = int(args.sort)

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
        sort_order = sorted_bin_list[:,-1]
        
        transition_matrix = sortMatrix(transition_matrix, sort_order)
        coordinate_ids_for_plot = sorted_bin_list[:,:-1]        
                
        
    # output
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


