#!/usr/bin/python2
import argparse
import sys
import numpy as np
from lib.logger import Logger
import lib.functions_ana_general as f
import lib.reweighting as reweighting
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def isStartState(_bin, start_state):
    if all(_bin.getCoordinateIds() == start_state):
        return True
    else:
        return False

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=
    'hdWE steady state analysis.')
parser.add_argument('-l', '--log', type=str, dest="logdir",
                    required=True, default="hdWE-log", metavar="DIR",
                    help="The logdir to load.")
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1, metavar='INT',
                    help="Last iteration to to read.")
parser.add_argument('-B', '--trace-flux-start', dest="first_ana_iteration",
                    type=int, default=0, metavar='INT', 
                    help='Iteration to start actual trace_flux analysis from.')
parser.add_argument('-s', '--start-state-bin-coordinate-id', dest="start_state",
                    required=True, type=int, nargs=3, 
                    help="Start state bin id for steady state flux calculation.")  
parser.add_argument('-t', '--target-state-bin-coordinate-id', dest="target_state",
                    required=True, type=int, nargs=3, 
                    help="Start state bin id for steady state flux calculation.")
parser.add_argument('-o', '--output', dest="output_file", 
                    required=False, type=str, default='ana_SS',
                    metavar = "FILE", 
                    help="output filename for A and B flux properties")
parser.add_argument('-r', '--reweighting-range', dest="reweighting_range",
                    type=float, default=0.5, metavar="FLOAT",
                    help="Fraction of previous iterations used for reweighting rate calculation.")  
parser.add_argument('-w', '--reweighting-iterations', dest="reweighting_iterations",
                    type=int, default=-1, metavar="INT",
                    help="Apply reweighting to first N of the read iterations.")
parser.add_argument('-N', '--pmf-bins', dest="pmf_bins",
                    type=int, default=100, metavar="INT",
                    help="Number of bins for the PMF.")
parser.add_argument('--tau', dest="tau",
                    type=float, default = 1.0, required=False,
                    help="Segment MD propagation time tau used in hdWE simulations.")
parser.add_argument('-p', '--plot', dest="plot", 
                    required=False, default=False, action="store_true",
                    help="Plot probabilities from start state per bin and iteration")
# Initialize
args = parser.parse_args()
first_iteration = 0
last_iteration  = args.last_iteration
logger = Logger(args.logdir)
if last_iteration < 0:
    last_iteration = logger.getLastIterationId()
# get the bin Ids of start state and target state
start_state  = args.start_state 
target_state = args.target_state
last_iteration_tmp = logger.loadIteration(last_iteration)
found_target_state = False
for this_bin in last_iteration_tmp.bins:
    if all(this_bin.getCoordinateIds() == start_state):
        start_state_id = this_bin.getId()
    if all(this_bin.getCoordinateIds() == target_state):
        target_state_id = this_bin.getId()
        found_target_state = True
if found_target_state == False:
    target_state_id = 1e99
    print('Warning: Target state not found in given iterations')
        
       
current_iteration = logger.loadIteration(first_iteration)

# assign initial probabilities
N = current_iteration.getNumberOfSegments()
for this_bin in current_iteration.bins:
    if isStartState(this_bin, start_state) == True:
        for this_segment in this_bin.segments:
            this_segment.setProbability(np.array([float(1.0/N),float(0.0)]))
        for this_segment in this_bin.initial_segments:
            this_segment.setProbability(np.array([float(1.0/N),float(0.0)]))           
    else:
        for this_segment in this_bin.segments:
            this_segment.setProbability(np.array([float(0.0),float(1.0/N)]))
        for this_segment in this_bin.initial_segments:
            this_segment.setProbability(np.array([float(0.0),float(1.0/N)]))  

flux = []
pmf_segment_data = []
from_start_state_probabilities = []

reweighter = reweighting.Reweighting(reweighting_range=args.reweighting_range)
reweighter.storeRateMatrix(current_iteration)

# Iteration Loop
for i in range(first_iteration + 1, last_iteration + 1):
    previous_iteration = current_iteration
    current_iteration = logger.loadIteration(i)
    sys.stdout.write('Iteration: {:05d}, Active Bins: {:05d}, Total Prob.: {}\r'.
                     format(i, current_iteration.getNumberOfActiveBins(), previous_iteration.getProbability()))
    sys.stdout.flush()

    # SPLITTING          
    for this_bin in current_iteration:
        # assign new probability to initial segments
        for this_initial_segment in this_bin.initial_segments:
            new_probability = previous_iteration\
                                .bins[this_initial_segment.getParentBinId()]\
                                .segments[this_initial_segment.getParentSegmentId()].getProbability()
            this_initial_segment.setProbability(new_probability)


        if this_bin.sample_region == True:
            if this_bin.getNumberOfSegments() == this_bin.getNumberOfInitialSegments():
                for this_segment in this_bin:
                    this_segment.setProbability(this_bin.initial_segments[this_segment.getId()].getProbability())
            elif this_bin.getNumberOfSegments() > this_bin.getNumberOfInitialSegments():
                split_dict = {}
                for this_segment in this_bin:
                    split_dict[this_segment.getParentNameString()] = split_dict.get(this_segment.getParentNameString() , 0) + 1
                
                for this_initial_segment in this_bin.initial_segments:
                    split_dict[this_initial_segment.getParentNameString()] = this_initial_segment.getProbability() / \
                                                    split_dict[this_initial_segment.getParentNameString()]
                for this_segment in this_bin:
                    this_segment.setProbability(split_dict[this_segment.getParentNameString()])
    
    # Apply the Steady State rules
    # Reaching the target state
    flux.append([])
    if current_iteration.getNumberOfBins() > target_state_id: 
        probability_from_start_state = 0.0
        for this_initial_segment in current_iteration.bins[target_state_id].initial_segments:
            # Only add probability to flux that comes from start state
            if this_initial_segment.getProbability()[0] > 0.0:
                flux[-1].append(this_initial_segment.getProbability()[0])
            # shift probability to target state history
            this_segment.setProbability( np.array([0, sum(this_segment.getProbability())] ) )
        
    # Reset Outer Region Bins
    current_iteration.resetOuterRegion(steady_state = True)

    # Reaching the start state
    for this_bin in current_iteration:
        if isStartState(this_bin, start_state) == True:
            for this_segment in this_bin:
                probability = this_segment.getProbability()
                this_segment.setProbability(np.array([sum(probability), 0.0]))
                
    # Set probability of start state to 1:
    start_state_probability = sum( current_iteration.bins[start_state_id].getProbability() )
    factor                  = 1.0 / start_state_probability
    for this_segment in current_iteration.bins[start_state_id]:
        this_segment.multProbability(factor)
    
    # Reweighting of bin probabilities
    if i < args.reweighting_iterations:
        # Keep track of the rate matrix
        reweighter.storeRateMatrix(current_iteration)
        if current_iteration.getNumberOfBins() > 1:
            reweighter.reweightBinProbabilities(current_iteration)
        factor = 1.0 / sum(current_iteration.bins[start_state_id].getProbability())
        for this_bin in current_iteration.bins:
            for this_segment in this_bin:
                this_segment.multProbability(factor)            
            
    # keep track of PMF-relevant segment data
    if i > args.first_ana_iteration:
        for this_bin in current_iteration:
            for this_segment in this_bin:
                #TODO: lazy 1d implementation
                parent_bin_id = this_segment.getParentBinId()
                parent_seg_id = this_segment.getParentSegmentId()  
                pmf_segment_data.append([ previous_iteration.bins[parent_bin_id].segments[parent_seg_id].getCoordinates()[0], 
                                          np.sum(this_segment.getProbability()) ])
    
    # Collect probability data for plotting
    if i > args.first_ana_iteration and args.plot: 
        iteration_bin_probs = [0.0] * (len(current_iteration.getBoundaries()[0]) + 1)
        for this_bin in current_iteration:
            if this_bin.getNumberOfInitialSegments() > 0:
                iteration_bin_probs[this_bin.getCoordinateIds()[0]] += this_bin.getInitialProbability()[0]
        iteration_bin_probs = np.apply_along_axis(np.log, 0, iteration_bin_probs)
        from_start_state_probabilities.append(iteration_bin_probs)

# Plot the probabilities from start state per bin
if args.plot:
    zdata = np.array(from_start_state_probabilities)
    xdata = range(zdata.shape[0])
    ydata = range(zdata.shape[1])
    x, y  = np.meshgrid(ydata, xdata)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel("Bin ID")
    ax.set_ylabel("Iteration")
    ax.set_zlabel("log(P(from start state))")
    ax.plot_wireframe(x, y, zdata, label="Probability from start state per bin")
    #ax.plot_surface(x, y, zdata, rstride=4, cstride=4, alpha=0.4)
    ax.legend()
    plt.show()

print '\n'
print 'flux: ', flux
for i in range(len(flux)):
    if flux[i] == []:
        flux[i] = [0.0]
    flux[i] = sum(flux[i])
print 'mean flux: {:8.7e}'.format( np.mean(flux[args.first_ana_iteration:]) / args.tau )

# print flux to file
fout = open(args.output_file + '.flux', 'w')
fout.write('# iteration   flux')
for i in range(len(flux)):
    fout.write('{:05d} {:8.7e}\n'.format(i, flux[i]))
fout.close()   

# calculate PMF
histogram = f.weightedHistogram(pmf_segment_data[args.first_ana_iteration:], args.pmf_bins)
fout = open(args.output_file + '.hist', 'w')
fout.write('# coord   probability')
for i in range(len(histogram)):
    fout.write('{:8.7e} {:8.7e}\n'.format(histogram[i,0], histogram[i,1]))
fout.close()
