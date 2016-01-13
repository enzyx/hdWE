#!/usr/bin/env python2
"""
Show the number of segments per iteration.
"""
from __future__ import print_function
import numpy
from lib.logger import Logger
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
parser.add_argument('-u', '--evolution-over-time', dest="evol", 
                    required=False, default=False, action="store_true",
                    help="Show the bin evolution over time (coordinate id 0).")
parser.add_argument('-p', '--plot', dest="plot", 
                    required=False, default=False, action="store_true",
                    help="Plot the fit with mathplotlib")


args = parser.parse_args()
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration, verbose=True)

########### Functions ##########
def showSegmentsPerIteration():
    xdata = []
    ydata = []
    
    for iteration in iterations:
        print("Iteration: {} #Segments: {: 6d}".format(iteration.getNameString(), iteration.getNumberOfSegments()))
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

def showBinEvolutionPerTime():
    zdata = []
    for iteration in iterations:
        segments_per_coordinate = [0] * (len(iteration.getBoundaries()[0]) + 1)
        for this_bin in iteration:
            segments_per_coordinate[this_bin.getCoordinateIds()[0]] += this_bin.getNumberOfSegments()
        if not args.plot:
            print(segments_per_coordinate)
        zdata.append(segments_per_coordinate)
    
    if args.plot:
        zdata = numpy.array(zdata)
        xdata = range(zdata.shape[0])
        ydata = range(zdata.shape[1])
        x, y  = numpy.meshgrid(ydata, xdata)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel("Bin ID")
        ax.set_ylabel("Iteration")
        ax.set_zlabel("#Segments")
        ax.plot_wireframe(x, y, zdata, label="#Segments/bin evolution")
        ax.legend()
        plt.show()


############# Main #############
if args.evol:
    showBinEvolutionPerTime()
else:
    showSegmentsPerIteration()