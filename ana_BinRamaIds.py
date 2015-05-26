#!/usr/bin/python2
import argparse
import sys
import ConfigParser
import matplotlib.pyplot as plt         ## plot library
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
import constants
from iteration import Iteration
from segment import Segment
from logger import Logger
import convergenceCheck
from hdWE_parameters import HdWEParameters


#### classes ####


            
##### functions #####


###### Parse command line ###### 
parser = argparse.ArgumentParser(description=
    'Bare Model to load iterations. ')
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0, metavar='INT',
                    help="First iteration to use.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1, metavar='INT',
                    help="Last iteration to to use.")  
parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    required=True, default="hdWE-log", metavar="FILE",
                    help="The logdir to load.")
parser.add_argument('-p', '--plot', dest="plot", 
                    default=False, action="store_true",
                    help="Plot the result.")
           
# Initialize
args = parser.parse_args()
                  
################################

# Get the iterations and parameters
logger = Logger(args.logdir)
iterations = logger.loadIterations(begin = args.first_iteration, 
                                   end   = args.last_iteration)


# a list of tuples (bin_id, rama_id, probability)
bin_datas = []    
for loop_bin in iterations[-1]:
    bin_datas.append( (loop_bin.getId(), loop_bin.getRamaId(), loop_bin.getProbability()) )

# sort according to probability
bin_datas = sorted(bin_datas, key=lambda x: x[2],reverse=True)
    
print ("Ramachandran Ids of bins and their probability:")
for bin_data in bin_datas:
    print ("{bin_id:05d}   {bin_rama_id}   {p}".format(
                bin_id      = bin_data[0],
                bin_rama_id = bin_data[1],
                p           = bin_data[2]))    
