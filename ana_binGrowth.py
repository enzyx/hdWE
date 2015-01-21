#!/usr/bin/python3
import argparse
import numpy
import constants
import sys
from logger import Logger
from math import log
from amber_module import MD_module

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=
    'Calculates the PMF along an arbitrary coordinate from the data of a hdWE run.')
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    required=False, type=int, default=0,
                    help="First iteration to use for PMF calculation.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    required=False, type=int, default=-1,
                    help="Last iteration to to use for PMF calculation.")  
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    required=True, default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")
################################

sys.stderr.write('\033[1mWriting numer of bins per iteration\033[0m\n')
               
# Initialize
args = parser.parse_args()

# Get the actual Iteration from logger module
logger = Logger(args.logfile)
iterations = logger.load_iterations(first = args.first_iteration, 
                                    last = args.last_iteration)
logger.close()

for iteration in iterations:
    print("{0: 5d} {1: 10d}".format(iteration.getId(),
                                    iteration.getNumberOfBins() ))
