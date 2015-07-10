#!/usr/bin/python2
import argparse
import sys
import numpy as np
from lib.logger import Logger
import lib.functions_ana_general as f
import lib.reweighting as reweighting
import lib.constants as constants
from math import log

#####################
###### classes ######
#####################

#####################
##### functions #####
#####################

def getStateFromCoordinate(segment, state_A, state_B):
    """
    Returns the state of a segment
    @return string state
    """
    state_per_dimension = []
    for coordinate in segment.getCoordinates():
        if coordinate > state_A[0] and coordinate <= state_A[1]:
            state_per_dimension.append('A')
        elif coordinate > state_B[0] and coordinate <= state_B[1]:
            state_per_dimension.append('B')
        else:
            state_per_dimension.append('0')
    
    this_state = state_per_dimension[0]
    for state in state_per_dimension[1:]:
            if state != this_state:
                this_state = '0'
                break 
    return this_state         

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=
    'Bare Model to load iterations. ')
parser.add_argument('-l', '--log', type=str, dest="logdir",
                    required=True, default="hdWE-log", metavar="DIR",
                    help="The logdir to load.")
parser.add_argument('-B', '--trace-flux-start', dest="trace_flux_start",
                    type=int, default=0, metavar='INT')
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0, metavar='INT',
                    help="First iteration to use.")
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1, metavar='INT',
                    help="Last iteration to to use.")
parser.add_argument('--state-A', dest="state_A",
                    required=False, type=float, nargs=2, 
                    help="Boundaries of the start state for rate calculation.")  
parser.add_argument('--state-B', dest="state_B",
                    required=False, type=float, nargs=2, 
                    help="Boundaries of the end state for rate calculation.")
parser.add_argument('-o', '--output', dest="output_file", 
                    required=False, type=str, default='ana_trace_flux',
                    help="output filename for A and B flux properties")
parser.add_argument('-r', '--reweighting-range', dest="reweighting_range",
                    type=float, default=0.5,
                    help="Reweighting range.")  
parser.add_argument('-w', '--reweighting-iterations', dest="reweighting_iterations",
                    type=int, default=-1,
                    help="Apply reweighting to first N iterations.")
parser.add_argument('-N', '--pmf-bins', dest="pmf_bins",
                    type=int, default=100,
                    help="Number of bins for the PMF..")

# Initialize
args = parser.parse_args()
first_iteration = args.trace_flux_start
last_iteration  = args.last_iteration
state_A         = np.sort(args.state_A) 
state_B         = np.sort(args.state_B)

################################

# Get the iterations and parameters
logger = Logger(args.logdir)
if last_iteration < 0:
    last_iteration = logger.getLastIterationId()

current_iteration = logger.loadIteration(first_iteration)
N = current_iteration.getNumberOfSegments()

# assign initial probabilities
for this_bin in current_iteration:
    for this_segment in this_bin:
        this_state = getStateFromCoordinate(this_segment, state_A, state_B)
        if this_state == 'A':
            this_segment.setProbability(np.array([1.0/N,0.0]))
        elif this_state == 'B':
            this_segment.setProbability(np.array([0.0,1.0/N]))
        else:
            this_segment.setProbability(np.array([0.5/N,0.5/N]))
    for this_segment in this_bin.initial_segments:
        if this_state == 'A':
            this_segment.setProbability(np.array([1.0/N,0.0]))
        elif this_state == 'B':
            this_segment.setProbability(np.array([0.0,1.0/N]))
        else:
            this_segment.setProbability(np.array([0.5/N,0.5/N]))
                     

flux_into_A         = []
flux_into_B         = []
probability_state_A = []
probability_state_B = []
probability_from_A  = []
probability_from_B  = []
pmf_segment_data    = []

reweighter = reweighting.Reweighting(reweighting_range=args.reweighting_range)
reweighter.storeRateMatrix(current_iteration)

# Iteration Loop
for i in range(first_iteration + 1, last_iteration + 1):
    sys.stdout.write('Iteration {:05d}\r'.format(i))
    sys.stdout.flush()
    previous_iteration = current_iteration
    current_iteration = logger.loadIteration(i)
    
    # Initialize data
    flux_into_A_iter         = 0.0 
    flux_into_B_iter         = 0.0
    probability_state_A_iter = 0.0 
    probability_state_B_iter = 0.0 
    probability_from_A_iter  = 0.0 
    probability_from_B_iter  = 0.0  
       
    for this_bin in current_iteration:
                    
        # assign new probability to initial segments
        for this_initial_segment in this_bin.initial_segments:
            new_probability = previous_iteration\
                                .bins[this_initial_segment.getParentBinId()]\
                                .segments[this_initial_segment.getParentSegmentId()].getProbability()
            this_initial_segment.setProbability(new_probability)

        # No split and merge
        if this_bin.getNumberOfSegments() == this_bin.getNumberOfInitialSegments():
            for this_segment in this_bin:
                this_segment.setProbability(this_bin.initial_segments[this_segment.getId()].getProbability())
        
        # Split
        elif this_bin.getNumberOfSegments() > this_bin.getNumberOfInitialSegments():
            split_dict = {}
            for this_segment in this_bin:
                try:
                    split_dict[this_segment.getParentNameString()] += 1
                except KeyError:
                    split_dict[this_segment.getParentNameString()] = 1
            
            for this_initial_segment in this_bin.initial_segments:
                split_dict[this_initial_segment.getParentNameString()] = this_initial_segment.getProbability() / \
                                                split_dict[this_initial_segment.getParentNameString()]
        
            for this_segment in this_bin:
                this_segment.setProbability(split_dict[this_segment.getParentNameString()])
        
        # Merge
        elif this_bin.getNumberOfSegments() < this_bin.getNumberOfInitialSegments():
            merge_dict = {}
            kept_prob = 0.0
            merged_probability = 0.0
            # Collect all initial probabilities in a dictionary
            for this_initial_segment in this_bin.initial_segments:
                merge_dict[this_initial_segment.getParentNameString()] = this_initial_segment.getProbability()
            # Assign initial probs to kept segments and set to 0.0 if consumed
            for this_segment in this_bin:
                this_segment.setProbability(merge_dict[this_segment.getParentNameString()])
                kept_prob += this_segment.getProbability()
                merge_dict[this_segment.getParentNameString()] = 0.0
            # Collect the remaining probability
            for parent_string_id in merge_dict:
                merged_probability += merge_dict[parent_string_id]

            merged_probability =  merged_probability/ this_bin.getNumberOfSegments()
            # assign merge probability

            for this_segment in this_bin:
                this_segment.setProbability( this_segment.getProbability() + merged_probability )
                
        # STATES
        for this_segment in this_bin:
            this_state = getStateFromCoordinate(this_segment, state_A, state_B)
            probability = this_segment.getProbability()
            if this_state == 'A':
                probability_state_A_iter += sum(probability)
                flux_into_A_iter  += probability[1]
                this_segment.setProbability(np.array([sum(probability), 0.0])) 
            elif this_state == 'B':
                probability_state_B_iter += sum(probability)
                flux_into_B_iter  += probability[0]
                this_segment.setProbability(np.array([0.0, sum(probability)]))   
                
    # keep track of PMF-relevant segment data
    if i > first_iteration:
        for this_bin in current_iteration:
            for this_segment in this_bin:
                #TODO: lazy 1d implementation
                pmf_segment_data.append([ this_segment.getCoordinates()[0], 
                                          np.sum(this_segment.getProbability()) ])
                     
    probability_from_A_iter = current_iteration.getProbability()[0]
    probability_from_B_iter = current_iteration.getProbability()[1]
                    
    flux_into_A.append(flux_into_A_iter)         
    flux_into_B.append(flux_into_B_iter) 
    probability_state_A.append(probability_state_A_iter)         
    probability_state_B.append(probability_state_B_iter)   
    probability_from_A.append(probability_from_A_iter)         
    probability_from_B.append(probability_from_B_iter)
    
    # Reweighting of bin probabilities
    #    The order of the following steps should no longer matter.  
    if i < args.reweighting_iterations:
        # Keep track of the rate matrix
        reweighter.storeRateMatrix(current_iteration)
        if current_iteration.getNumberOfBins() > 1:
            reweighter.reweightBinProbabilities(current_iteration)


##########################
######### OUTPUT #########
##########################
b = args.first_iteration - args.trace_flux_start
e = args.last_iteration

# Data time series
fout = open(args.output_file + '.dat', 'w')
fout.write("# F->A           P_A            P->A            F->B            P_B             P->B \n")
for i in range(len(flux_into_A)):
    fout.write("{:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}\n".format(
                         flux_into_A[i], 
                         probability_state_A[i], 
                         probability_from_A[i], 
                         flux_into_B[i], 
                         probability_state_B[i], 
                         probability_from_B[i]))
fout.close()

# Cumulative mean of the time series
cum_flux_into_A         = f.cumulative_mean(flux_into_A[b:e])
cum_probability_state_A = f.cumulative_mean(probability_state_A[b:e])
cum_probability_from_A  = f.cumulative_mean(probability_from_A[b:e])
cum_flux_into_B         = f.cumulative_mean(flux_into_B[b:e])
cum_probability_state_B = f.cumulative_mean(probability_state_B[b:e])
cum_probability_from_B  = f.cumulative_mean(probability_from_B[b:e])

fout = open(args.output_file + '.cum', 'w')
for i in range(len(cum_flux_into_A)):
    fout.write("{:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}\n".format(
                         cum_flux_into_A[i], 
                         cum_probability_state_A[i], 
                         cum_probability_from_A[i], 
                         cum_flux_into_B[i], 
                         cum_probability_state_B[i], 
                         cum_probability_from_B[i]))
fout.close()

# Data autocorrelaction functions
auto_flux_into_A         = f.autocorrelation_function(flux_into_A[b:e])
auto_probability_state_A = f.autocorrelation_function(probability_state_A[b:e])
auto_probability_from_A  = f.autocorrelation_function(probability_from_A[b:e])
auto_flux_into_B         = f.autocorrelation_function(flux_into_B[b:e])
auto_probability_state_B = f.autocorrelation_function(probability_state_B[b:e])
auto_probability_from_B  = f.autocorrelation_function(probability_from_B[b:e])

fout = open(args.output_file + '.auto', 'w')
for i in range(len(auto_flux_into_A)):
    fout.write("{:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}\n".format(
                         auto_flux_into_A[i], 
                         auto_probability_state_A[i], 
                         auto_probability_from_A[i], 
                         auto_flux_into_B[i], 
                         auto_probability_state_B[i], 
                         auto_probability_from_B[i]))
fout.close()

# calculate PMF
histogram = f.weightedHistogram(pmf_segment_data, args.pmf_bins)
pmf = np.zeros(len(histogram), dtype=float)
for i in range(len(histogram)):
    if histogram[i,1] > 0.0:
        pmf[i] = -constants.kT *  log(histogram[i,1])
pmf -= np.min(pmf)

fout = open(args.output_file + '.pmf', 'w')
fout.write('# coord   free energy   probability')
for i in range(len(pmf)):
    fout.write('{:8.7e} {:8.7e} {:8.7e}\n'.format(histogram[i,0], pmf[i], histogram[i,1]))
fout.close()

# Output
block_size = 10

print "A --> B"
print   f.block_bootstrap(flux_into_B[b:e], np.mean, block_size)    
print   f.block_bootstrap(probability_state_A[b:e], np.mean, block_size)
print   f.block_bootstrap(probability_from_A[b:e], np.mean, block_size)  
print 'rate:'
print   np.mean(flux_into_B[b:e])  / np.mean(probability_state_A[b:e])
print '1/mfpt:'
print   np.mean(flux_into_B[b:e])  / np.mean(probability_from_A[b:e])

print "B --> A"
print   f.block_bootstrap(flux_into_A[b:e], np.mean, block_size)    
print   f.block_bootstrap(probability_state_B[b:e], np.mean, block_size)
print   f.block_bootstrap(probability_from_B[b:e], np.mean, block_size)  
print 'rate:'
print   np.mean(flux_into_A[b:e])  / np.mean(probability_state_B[b:e])
print '1/mfpt:'
print   np.mean(flux_into_A[b:e])  / np.mean(probability_from_B[b:e])    
