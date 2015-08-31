#!/usr/bin/python2
import argparse
import sys
import numpy as np
from lib.logger import Logger
import lib.functions_ana_general as f
import lib.reweighting as reweighting
import lib.constants as constants
from math import log

def isStartState(bin, start_state):
    if all(bin.getCoordinateIds() == start_state):
        return True
    else:
        return False

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=
    'hdWE steady state analysis.')
parser.add_argument('-l', '--log', type=str, dest="logdir",
                    required=True, default="hdWE-log", metavar="DIR",
                    help="The logdir to load.")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0, metavar='INT',
                    help="First iteration to read.")
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1, metavar='INT',
                    help="Last iteration to to read.")
parser.add_argument('-B', '--trace-flux-start', dest="first_ana_iteration",
                    type=int, default=0, metavar='INT', 
                    help='Iteration to start actual trace_flux analysis from.')
parser.add_argument('-start', '--start-state-bin-coordinate-id', dest="start_state",
                    required=True, type=int, nargs=3, 
                    help="Start state bin id for steady state flux calculation.")  
parser.add_argument('-target', '--target-state-bin-coordinate-id', dest="target_state",
                    required=True, type=int, nargs=3, 
                    help="Start state bin id for steady state flux calculation.")
parser.add_argument('-o', '--output', dest="output_file", 
                    required=False, type=str, default='ana_trace_flux',
                    metavar = "FILE", 
                    help="output filename for A and B flux properties")
parser.add_argument('-r', '--reweighting-range', dest="reweighting_range",
                    type=float, default=0.5, metavar="FLOAT",
                    help="Fraction of previous iterations used for reweighting rate calculation.")  
parser.add_argument('-w', '--reweighting-iterations', dest="reweighting_iterations",
                    type=int, default=-1, metavar="INT",
                    help="Apply reweighting to first N of the read iterations.")
parser.add_argument('-N', '--pmf-bins', dest="pmf_bins",
                    type=int, default=200, metavar="INT",
                    help="Number of bins for the PMF.")
parser.add_argument('-t', '--tau', dest="tau",
                    type=float, default = 0, required=False,
                    help="tau used in hdWE simulations.")

# Initialize
args = parser.parse_args()
first_iteration = args.first_iteration
last_iteration  = args.last_iteration
logger = Logger(args.logdir)
if last_iteration < 0:
    last_iteration = logger.getLastIterationId()
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
    print('Target state not found within given range of iterations')
    sys.exit()
        
        
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

reweighter = reweighting.Reweighting(reweighting_range=args.reweighting_range)
reweighter.storeRateMatrix(current_iteration)

# Iteration Loop
for i in range(first_iteration + 1, last_iteration + 1):
    previous_iteration = current_iteration
    current_iteration = logger.loadIteration(i)
    sys.stdout.write('Iteration: {:05d}, Active Bins: {:05d}, Total Prob.: {:1.8f}\r'.
                     format(i, current_iteration.getNumberOfActiveBins(), sum(previous_iteration.getProbability())))
    sys.stdout.flush()
          
    for this_bin in current_iteration:
                    
        # assign new probability to initial segments
        for this_initial_segment in this_bin.initial_segments:
            new_probability = previous_iteration\
                                .bins[this_initial_segment.getParentBinId()]\
                                .segments[this_initial_segment.getParentSegmentId()].getProbability()
            this_initial_segment.setProbability(new_probability)


        if this_bin.sample_region == True:
            # No Split
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

    # Reset Outer Region Bins
    current_iteration.resetOuterRegion(target_state_id = target_state_id)
    
    # Reaching the start state
    for this_bin in current_iteration:
        if isStartState(this_bin, start_state) == True:
            for this_segment in this_bin:
                probability = this_segment.getProbability()
                this_segment.setProbability(np.array([sum(probability), 0.0]))
               
    
    # Incoming flux
    if current_iteration.getNumberOfBins() > target_state_id: 
        probability_from_start_state = 0.0
        # empty target bin has no tuple probability
        if all(current_iteration.bins[target_state_id].getInitialProbability())  == 0.0:
            total_probability = np.array([0.0, 0.0])
        else:
            total_probability = sum( current_iteration.bins[target_state_id].getInitialProbability() )
        for this_segment in current_iteration.bins[target_state_id].segments:
            # Only add probability that comes from start state
            probability_from_start_state += this_segment.getProbability()[0]
        flux.append(probability_from_start_state)
        # delete the segments
        current_iteration.bins[target_state_id].deleteAllSegments()
        # recycle the total probability
        N_target_segments = current_iteration.bins[target_state_id].getNumberOfSegments()
        reclycled_probability_per_segment = 1.0 * N_target_segments / total_probability
        for this_segment in current_iteration.bins[target_state_id].segments:
            this_segment.addProbability([reclycled_probability_per_segment, 0.0])
        
        
    
    
    # Reweighting of bin probabilities
    if i < args.reweighting_iterations:
        # Keep track of the rate matrix
        reweighter.storeRateMatrix(current_iteration)
        if current_iteration.getNumberOfBins() > 1:
            reweighter.reweightBinProbabilities(current_iteration)

print '\nflux: ', flux
    

##########################
######### OUTPUT #########
##########################
b = args.first_ana_iteration - args.first_iteration
e = args.last_iteration

# # Data time series
# fout = open(args.output_file + '.dat', 'w')
# fout.write("# F->A           P_A            P->A            F->B            P_B             P->B \n")
# for i in range(len(flux_into_A)):
#     fout.write("{:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}\n".format(
#                          flux_into_A[i], 
#                          probability_state_A[i], 
#                          probability_from_A[i], 
#                          flux_into_B[i], 
#                          probability_state_B[i], 
#                          probability_from_B[i]))
# fout.close()
# 
# # Cumulative mean of the time series
# cum_flux_into_A         = f.cumulative_mean(flux_into_A[b:e])
# cum_probability_state_A = f.cumulative_mean(probability_state_A[b:e])
# cum_probability_from_A  = f.cumulative_mean(probability_from_A[b:e])
# cum_flux_into_B         = f.cumulative_mean(flux_into_B[b:e])
# cum_probability_state_B = f.cumulative_mean(probability_state_B[b:e])
# cum_probability_from_B  = f.cumulative_mean(probability_from_B[b:e])
# 
# fout = open(args.output_file + '.cum', 'w')
# for i in range(len(cum_flux_into_A)):
#     fout.write("{:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}\n".format(
#                          cum_flux_into_A[i], 
#                          cum_probability_state_A[i], 
#                          cum_probability_from_A[i], 
#                          cum_flux_into_B[i], 
#                          cum_probability_state_B[i], 
#                          cum_probability_from_B[i]))
# fout.close()
# 
# # Data autocorrelaction functions
# auto_flux_into_A         = f.autocorrelation_function(flux_into_A[b:e])
# auto_probability_state_A = f.autocorrelation_function(probability_state_A[b:e])
# auto_probability_from_A  = f.autocorrelation_function(probability_from_A[b:e])
# auto_flux_into_B         = f.autocorrelation_function(flux_into_B[b:e])
# auto_probability_state_B = f.autocorrelation_function(probability_state_B[b:e])
# auto_probability_from_B  = f.autocorrelation_function(probability_from_B[b:e])
# 
# fout = open(args.output_file + '.auto', 'w')
# for i in range(len(auto_flux_into_A)):
#     fout.write("{:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}\n".format(
#                          auto_flux_into_A[i], 
#                          auto_probability_state_A[i], 
#                          auto_probability_from_A[i], 
#                          auto_flux_into_B[i], 
#                          auto_probability_state_B[i], 
#                          auto_probability_from_B[i]))
# fout.close()
# 
# # calculate PMF
# histogram = f.weightedHistogram(pmf_segment_data, args.pmf_bins)
# pmf = np.zeros(len(histogram), dtype=float)
# for i in range(len(histogram)):
#     if histogram[i,1] > 0.0:
#         pmf[i] = -constants.kT *  log(histogram[i,1])
#     else:
#         pmf[i] = 'Inf'
# pmf -= np.min(pmf)
# 
# fout = open(args.output_file + '.pmf', 'w')
# fout.write('# coord   free energy   probability')
# for i in range(len(pmf)):
#     fout.write('{:8.7e} {:8.7e} {:8.7e}\n'.format(histogram[i,0], pmf[i], histogram[i,1]))
# fout.close()
# 
# # Output
# block_size = 10
# 
# ############################
# #      Pretty print        #
# ############################
# print("States:  A=[{},{}]  B=[{},{}]".format(state_A[0], state_A[1], state_B[0], state_B[1]))
# # State A 
# print
# print(" State A -> B:")
# block_bootstrap_flux_into_B =  f.block_bootstrap(flux_into_B[b:e], np.mean, block_size)
# print("    Flux(A->B):  {:5.4e}  CI: [{:5.4e}, {:5.4e}]".format(block_bootstrap_flux_into_B[0], 
#                                                                block_bootstrap_flux_into_B[1][0],
#                                                                block_bootstrap_flux_into_B[1][1]))
# block_bootstrap_prop_A = f.block_bootstrap(probability_state_A[b:e], np.mean, block_size)
# print("    P(A):        {:5.4e}  CI: [{:5.4e}, {:5.4e}]".format(block_bootstrap_prop_A[0],
#                                                                block_bootstrap_prop_A[1][0],
#                                                                block_bootstrap_prop_A[1][1]))
# block_bootstrap_prop_from_A = f.block_bootstrap(probability_from_A[b:e], np.mean, block_size)
# print("    P(from A):   {:5.4e}  CI: [{:5.4e}, {:5.4e}]".format(block_bootstrap_prop_from_A[0],
#                                                                block_bootstrap_prop_from_A[1][0],
#                                                                block_bootstrap_prop_from_A[1][1]))
# if args.tau > 0:
#     print("    k(A->B):     {:5.4e} 1/s".format((np.mean(flux_into_B[b:e])  / np.mean(probability_state_A[b:e]))/args.tau))
# else:
#     print("    k(A->B):     {:5.4e}".format(np.mean(flux_into_B[b:e])  / np.mean(probability_state_A[b:e])))
# print("    1/(MFPT):    {:5.4e}".format(np.mean(flux_into_B[b:e])  / np.mean(probability_from_A[b:e])))
# 
# # State B
# print
# print(" State B -> A:")
# block_bootstrap_flux_into_A =  f.block_bootstrap(flux_into_A[b:e], np.mean, block_size)
# print("    Flux(B->A):  {:5.4e}  CI: [{:5.4e}, {:5.4e}]".format(block_bootstrap_flux_into_A[0], 
#                                                                block_bootstrap_flux_into_A[1][0],
#                                                                block_bootstrap_flux_into_A[1][1]))
# block_bootstrap_prop_B = f.block_bootstrap(probability_state_B[b:e], np.mean, block_size)
# print("    P(B):        {:5.4e}  CI: [{:5.4e}, {:5.4e}]".format(block_bootstrap_prop_B[0],
#                                                                block_bootstrap_prop_B[1][0],
#                                                                block_bootstrap_prop_B[1][1]))
# block_bootstrap_prop_from_B = f.block_bootstrap(probability_from_B[b:e], np.mean, block_size)
# print("    P(from B):   {:5.4e}  CI: [{:5.4e}, {:5.4e}]".format(block_bootstrap_prop_from_B[0],
#                                                                block_bootstrap_prop_from_B[1][0],
#                                                                block_bootstrap_prop_from_B[1][1]))
# if args.tau > 0:
#     print("    k(B->A):     {:5.4e} 1/s".format((np.mean(flux_into_A[b:e])  / np.mean(probability_state_B[b:e]))/args.tau))
# else:
#     print("    k(B->A):     {:5.4e}".format(np.mean(flux_into_A[b:e])  / np.mean(probability_state_B[b:e])))
# print("    1/(MFPT):    {:5.4e}".format(np.mean(flux_into_A[b:e])  / np.mean(probability_from_B[b:e])))
# 
# ########################################

