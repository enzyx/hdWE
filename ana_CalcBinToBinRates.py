#!/usr/bin/env python2
"""

"""


from __future__ import print_function
import numpy
from logger import Logger
import argparse  

###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")                  
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0,
                    help="First iteration to use for PMF calculation.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Last iteration to to use for PMF calculation.")  
parser.add_argument('-t', '--bin-to-bin-transitions', action="store_true",
                    dest="bin_to_bin_transitions", 
                    help="Print a matrix with bin to bin transition events")


args = parser.parse_args()
#get the actual Iteration from logger module
logger = Logger(args.logfile, APPEND = False)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration)
logger.close()

for iteration in iterations:
    print("Iteration: {0: 5d}".format(iteration.getId()))
    for _bin in iterations[-1]:
        print("Bin: {0: 5d} p={p:g}".format(_bin.getId(), p=_bin.getProbability()))
        for segment in _bin:
            print("    Segment {0: 5d} p = {p:g}: parent segment {1: 5d} bin {2: 5d}".format(segment.getId(),
                                            segment.getParentSegmentId(),
                                            segment.getParentBinId(),
                                            p=segment.getProbability()))

