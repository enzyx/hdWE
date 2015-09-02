#!/usr/bin/env python2
"""
Show the number of segments per iteration.
"""
from __future__ import print_function
import numpy
from lib.logger import Logger
import argparse
import matplotlib.pyplot as plt

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
parser.add_argument('-p', '--plot', dest="plot", 
                    required=False, default=False, action="store_true",
                    help="Plot the fit with mathplotlib")

args = parser.parse_args()
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration, verbose=True)

xdata = []
ydata = []

for iteration in iterations:
    print("Iteration: {:05d} #Segments: {: 6d}".format(iteration.getId(), iteration.getNumberOfSegments()))
    xdata.append(iteration.getId())
    ydata.append(iteration.getNumberOfSegments())

# Plot if required
if args.plot:
    xdata = numpy.array(xdata)
    ydata = numpy.array(ydata)
    plt.plot(xdata, ydata, label="#Segments per Iteration")
    plt.ylabel("#Segments")
    plt.xlabel("Iteration")
    plt.legend()
    plt.show()

