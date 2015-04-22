#!/usr/bin/env python2
"""

"""


from __future__ import print_function
import numpy
from logger import Logger
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
parser.add_argument('-t', '--bin-to-bin-transitions', action="store_true",
                    dest="bin_to_bin_transitions", 
                    help="Print a matrix with bin to bin transition events")


args = parser.parse_args()
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration)

N = iterations[-1].getNumberOfBins()
transition_matrix = numpy.zeros((N,N))

for iteration in iterations:
    for this_bin in iteration:
        for segment in this_bin.initial_segments:
            transition_matrix[segment.getParentBinId(), segment.getBinId()] += 1

print( transition_matrix )

