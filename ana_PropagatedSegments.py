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
parser.add_argument('-v', '--verbose', action="store_true",
                    dest="verbose", 
                    help="Print information for every iteration")

args = parser.parse_args()
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration, verbose=True)

VERBOSE         = args.verbose  

total_count = 0
for iteration in iterations:
    total_iteration_segments = 0
    active_bins = 0
    inactive_bins = 0
    for this_bin in iteration:
        this_number_of_segments = this_bin.getNumberOfSegments()
        if this_bin.getSampleRegion():
            active_bins += 1
            total_count += this_number_of_segments
            total_iteration_segments += this_number_of_segments
        else:
            inactive_bins += 1
    if VERBOSE:
        print("Iteration: {:> 5d} active (inactive) #bins: {:> 5d} ({:> 5d}) #segments: {:> 5d}".format(iteration.getId(), 
                                                                         active_bins,
                                                                         inactive_bins,
                                                                         total_iteration_segments))

print("""---------------------------------------
Total propagated #segments: {:> 8d}
---------------------------------------""".format(total_count))
