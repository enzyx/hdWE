#!/usr/bin/python2
from __future__ import print_function
from logger import Logger
import argparse 

parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Last iteration to use.") 
                    
# Initialize
print('\033[1mSumming up all propagated segments (without closed bins).\033[0m')      
args = parser.parse_args()


#get the actual Iteration from logger module
logger = Logger(args.logfile, append = True)
iterations = logger.loadIterations(0, args.last_iteration)
logger.close()

n_segments = 0
for iteration_loop in iterations:
    n_segments += iteration_loop.getNumberOfPropagatedSegments()
print('Total number of propagated segments: ' + str(n_segments-1))