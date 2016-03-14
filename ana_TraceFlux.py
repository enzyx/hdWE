#!/usr/bin/python2
import argparse
import sys
import numpy as np
from lib.logger import Logger
import lib.functions_ana_general as f
import lib.reweighting as reweighting
import lib.constants as constants
import math
import matplotlib.pyplot as plt
import scipy.integrate

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=
    'Bare Model to load iterations. ')
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
parser.add_argument('--state-A', dest="state_A", metavar="FLOAT",
                    required=False, type=float, nargs=2, 
                    help="Boundaries of the start state for rate calculation.")  
parser.add_argument('--state-B', dest="state_B", metavar="FLOAT",
                    required=False, type=float, nargs=2, 
                    help="Boundaries of the end state for rate calculation.")
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
parser.add_argument('--velocities', dest="velocities",
                    action='store_true', required=False,
                    help="check velocity distributions and print them.")
parser.add_argument('--auto', dest="auto", default=False, 
                    const=1, nargs='?', required=False, type=int,
                    help="save autocorrelaction function. "\
                         "given value is frequency of data usage to speed up N^2 calculations")
parser.add_argument('-k', dest="k",
                    required=False, default = 0,type=float,
                    help="Force constant for analytical calculation of probabilities of the states for the two particle test system.")
parser.add_argument('-bs', dest="bs",
                    required=True,type=int,
                    help="bootstrap block size.")
parser.add_argument('--bs-samples', dest="bs_samples",
                    required=False,type=int, default=10000,
                    help="number of bootstrap samples.")
parser.add_argument('--rates-only', dest="rates_only",
                    action='store_true', required=False,
                    help="only calculate rates for output")


# Initialize
args = parser.parse_args()
first_iteration = args.first_iteration
last_iteration  = args.last_iteration
logger = Logger(args.logdir)
if last_iteration < 0:
    last_iteration = logger.getLastIterationId()

current_iteration = logger.loadIteration(first_iteration)

state_A         = np.sort(args.state_A) 
state_B         = np.sort(args.state_B)


# assign initial probabilities
N = current_iteration.getNumberOfSegments()
for this_bin in current_iteration:
    for this_segment in this_bin:
        this_state = f.getStateFromCoordinate(this_segment, state_A, state_B)
        if this_state == 'A':
            this_segment.setProbability(np.array([float(1.0/N),float(0.0),float(0.0)]))
        elif this_state == 'B':
            this_segment.setProbability(np.array([float(0.0),float(1.0/N),float(0.0)]))
        else:
            this_segment.setProbability(np.array([float(0.0),float(0.0),float(1.0/N)]))
    for this_segment in this_bin.initial_segments:
        if this_state == 'A':
            this_segment.setProbability(np.array([float(1.0/N),float(0.0),float(0.0)]))
        elif this_state == 'B':
            this_segment.setProbability(np.array([float(0.0),float(1.0/N),float(0.0)]))
        else:
            this_segment.setProbability(np.array([float(0.0),float(0.0),float(1.0/N)]))
                     

flux_into_A         = []
flux_into_B         = []
probability_state_A = []
probability_state_B = []
probability_from_A  = []
probability_from_B  = []
pmf_segment_data    = []
pmf_segment_data_A  = []
pmf_segment_data_B  = []
velocity_data          = {'A': [], 'B':[]}
"""
stores all velocity_data from A and from B:
velocity_data['A'] is a list of datapoints [prob, coordinate, velocity]
"""

reweighter = reweighting.Reweighting(reweighting_range=args.reweighting_range)
reweighter.storeRateMatrix(current_iteration)

bin_prob_out = open('ana_trace_flux.binprobs.dat', 'w')
# Iteration Loop
merge_counter = 0
for i in range(first_iteration + 1, last_iteration + 1):
    previous_iteration = current_iteration
    current_iteration = logger.loadIteration(i)
    sys.stderr.write('Iteration: {:08d}, Active Bins: {:05d}, Total Prob.: {:1.8f}\r'.
                     format(i, current_iteration.getNumberOfActiveBins(), sum(previous_iteration.getProbability())))
    sys.stderr.flush()
        
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

        if this_bin.sample_region == True:
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
                # Skip merging if no merge_list exists, which means merge_mode was none
                if len(this_bin.merge_list) >= 1:
                    probabilities = []
                    for this_initial_segment in this_bin.initial_segments:
                        probabilities.append(this_initial_segment.getProbability())
                    
                    # strange behavior with list?  
                    probabilities = np.array(probabilities)
                    # reconstruct probabilities from merge list
                    for merge_entry in this_bin.merge_list:
                        # get merged probability
                        merged_probability = probabilities[merge_entry[0]]
                        merged_probability  = 1.0 * merged_probability / (len(merge_entry) - 1 )
                        for target_segment in merge_entry[1:]:
                            probabilities[target_segment] += np.array(merged_probability)                 
                        probabilities = np.delete(probabilities, merge_entry[0],0)
                    
                    # set merged probabilities of the segments 
                    for segment_index in range(0,len(this_bin.segments)):
                        this_bin.segments[segment_index].setProbability(probabilities[segment_index])
            
    # Reset Outer Region Bins
    current_iteration.resetOuterRegion(steady_state=True)

    # STATES
    for this_bin in current_iteration:
        for this_segment in this_bin:
            this_state = f.getStateFromCoordinate(this_segment, state_A, state_B)
            probability = this_segment.getProbability()
            if this_state == 'A':
                probability_state_A_iter += sum(probability)
                flux_into_A_iter  += probability[1]
                this_segment.setProbability(np.array([sum(probability), 0.0, 0.0])) 
            elif this_state == 'B':
                probability_state_B_iter += sum(probability)
                flux_into_B_iter  += probability[0]
                this_segment.setProbability(np.array([0.0, sum(probability), 0.0]))   

    # Reweighting of bin probabilities
    #    The order of the following steps should no longer matter.  
    if i < args.reweighting_iterations:
        # Keep track of the rate matrix
        reweighter.storeRateMatrix(current_iteration)
        if current_iteration.getNumberOfBins() > 1:
            reweighter.reweightBinProbabilities(current_iteration)
    
    # keep track of PMF-relevant segment data
    if i > args.first_ana_iteration:
        for this_bin in current_iteration:
            for this_segment in this_bin:
                #TODO: lazy 1d implementation
                #pmf_segment_data.append([ this_segment.getCoordinates()[0], 
                #                          np.sum(this_segment.getProbability()) ])
                # the before implementation took probabilities that corresponded to wrong coordinate values
                parent_bin_id = this_segment.getParentBinId()
                parent_seg_id = this_segment.getParentSegmentId()
                pmf_segment_data.append([ previous_iteration.bins[parent_bin_id].segments[parent_seg_id].getCoordinates()[0], 
                                          np.sum(this_segment.getProbability()) ])
                # PMF from A, from B
                pmf_segment_data_A.append([ previous_iteration.bins[parent_bin_id].segments[parent_seg_id].getCoordinates()[0], 
                                          this_segment.getProbability()[0] ])
                pmf_segment_data_B.append([ previous_iteration.bins[parent_bin_id].segments[parent_seg_id].getCoordinates()[0], 
                                          this_segment.getProbability()[1] ])
    
    probability_from_A_iter = current_iteration.getProbability()[0]
    probability_from_B_iter = current_iteration.getProbability()[1]
                    
    flux_into_A.append(flux_into_A_iter)         
    flux_into_B.append(flux_into_B_iter) 
    probability_state_A.append(probability_state_A_iter)         
    probability_state_B.append(probability_state_B_iter)   
    probability_from_A.append(probability_from_A_iter)         
    probability_from_B.append(probability_from_B_iter)
    
    
    for this_bin in current_iteration:
        prob = this_bin.getProbability()
        if type(prob) != float:
            prob = sum(prob)
        bin_prob_out.write("{: 8.7e}".format(prob))
    bin_prob_out.write('\n')
    bin_prob_out.flush()
    
    # store velocities
    if args.velocities:
        for this_bin in current_iteration:
            for this_segment in this_bin:      
                velocity_data['A'].append([this_segment.getProbability()[0], 
                                        this_segment.getCoordinates()[0], 
                                        this_segment.getVelocities()[0]])
                velocity_data['B'].append([this_segment.getProbability()[1], 
                                        this_segment.getCoordinates()[0], 
                                        this_segment.getVelocities()[0]])

bin_prob_out.close()

##########################
######### OUTPUT #########
##########################
sys.stderr.write('\n')
b = args.first_ana_iteration - args.first_iteration
e = last_iteration

# Data time series
sys.stderr.write('- writing time series to .dat\n')
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
sys.stderr.write('- calculating cumulative means\n')
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

# Data autocorrelation functions
if args.auto:
    a_every = args.auto
    sys.stderr.write('- calculating autocorrelations\n')
    auto_flux_into_A         = f.autocorrelation_function(flux_into_A[b:e:a_every])
    auto_probability_state_A = f.autocorrelation_function(probability_state_A[b:e:a_every])
    auto_probability_from_A  = f.autocorrelation_function(probability_from_A[b:e:a_every])
    auto_flux_into_B         = f.autocorrelation_function(flux_into_B[b:e:a_every])
    auto_probability_state_B = f.autocorrelation_function(probability_state_B[b:e:a_every])
    auto_probability_from_B  = f.autocorrelation_function(probability_from_B[b:e:a_every])
    
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
sys.stderr.write('- calculating PMFs\n')
histogram = f.weightedHistogram(pmf_segment_data, args.pmf_bins)
pmf = np.zeros(len(histogram), dtype=float)
for i in range(len(histogram)):
    if histogram[i,1] > 0.0:
        pmf[i] = -constants.kT *  math.log(histogram[i,1])
    else:
        pmf[i] = 'Inf'
pmf -= np.min(pmf)
  
fout = open(args.output_file + '.pmf', 'w')
fout.write('# coord   free energy   probability')
for i in range(len(pmf)):
    fout.write('{:8.7e} {:8.7e} {:8.7e}\n'.format(histogram[i,0], pmf[i], histogram[i,1]))
fout.close()

# calculate PMF from A 
histogram = f.weightedHistogram(pmf_segment_data_A, args.pmf_bins)
pmf = np.zeros(len(histogram), dtype=float)
for i in range(len(histogram)):
    if histogram[i,1] > 0.0:
        pmf[i] = -constants.kT *  math.log(histogram[i,1])
    else:
        pmf[i] = 'Inf'
pmf -= np.min(pmf)

fout = open(args.output_file + '_A.pmf', 'w')
fout.write('# coord   free energy   probability')
for i in range(len(pmf)):
    fout.write('{:8.7e} {:8.7e} {:8.7e}\n'.format(histogram[i,0], pmf[i], histogram[i,1]))
fout.close()

# calculate PMF from B
histogram = f.weightedHistogram(pmf_segment_data_B, args.pmf_bins)
pmf = np.zeros(len(histogram), dtype=float)
for i in range(len(histogram)):
    if histogram[i,1] > 0.0:
        pmf[i] = -constants.kT *  math.log(histogram[i,1])
    else:
        pmf[i] = 'Inf'
pmf -= np.min(pmf)

fout = open(args.output_file + '_B.pmf', 'w')
fout.write('# coord   free energy   probability')
for i in range(len(pmf)):
    fout.write('{:8.7e} {:8.7e} {:8.7e}\n'.format(histogram[i,0], pmf[i], histogram[i,1]))
fout.close()

#### calculate velocity histograms ####
if args.velocities:
    sys.stderr.write('- calculating velocity histograms\n')    
    # save velocity data to file
#     for state in ['A', 'B']:
#         v_out = open('ana_trace_flux.velo.{}.dat'.format(state), 'w')
#         v_out.write('#prob coord velocity\n')
#         for velo_entry in velocity_data[state]:
#             v_out.write('{} {} {}\n'.format(velo_entry[0], velo_entry[1], velo_entry[2]))
#         v_out.close()
    # print binned velocity histograms
    # plotting parameters
    colors={'A': 'green', 'B': 'blue'}
    n_histo_bins = 100
    v_range = 30
    histograms = []
    histograms = {'A': [], 'B': []}
    
    # bin-sorting paramters
    bin_boundaries = current_iteration.getBoundaries()[0]
    N_BINS = len(bin_boundaries) + 1
    binned_velocity_data = {}
    for state in ['A', 'B']:
        binned_velocity_data.update({state: []})
        for i in range(N_BINS):
            binned_velocity_data[state].append([])
        velocity_data[state].sort(key=lambda x: x[1])
        for velo_entry in velocity_data[state]:
            coord = velo_entry[1]
            for cbin_id, bin_boundary in enumerate(bin_boundaries):
                if coord < bin_boundary:
                    binned_velocity_data[state][cbin_id].append([velo_entry[2], velo_entry[0]]) # append velocity and probabililty
                    break
            # if it was higher than any boundary put into last bin
            if coord > bin_boundaries[-1]:
                binned_velocity_data[state][-1].append([velo_entry[2], velo_entry[0]])

    # create histograms
    for cbin_id in range(N_BINS):
        for state in ['A', 'B']:
            velo_data = np.asarray(binned_velocity_data[state][cbin_id])
            # foolproof for empty bins
            if len(velo_data) == 0:
                velo_data = np.asarray([[0.0, 0.0]])
            hist, bins = np.histogram(velo_data[:,0], 
                                      bins=n_histo_bins, 
                                      range = (-v_range, v_range),
                                      weights = velo_data[:,1],
                                      density = True)
            width = 1.0 * (bins[1] - bins[0])
            centers = (bins[:-1] + bins[1:]) / 2
            histograms[state].append((hist, centers))
    
    # save histograms
    hist_to_save = [np.asarray(histograms['A'][0][1])]
    for i in range(N_BINS):
        for state in ['A', 'B']:
            hist_to_save.append(np.asarray(histograms[state][i][0]))
    np.savetxt('ana_trace_flux.velo.histograms.dat', np.transpose(hist_to_save))

    # plotting
    for cbin_id in range(N_BINS):
        plt.clf()
        for state in ['A', 'B']:
            hist, centers = histograms[state][cbin_id]
#             n2, bins2, patches2 = plt.hist(velo_data[:,0], 
#                                            n_histo_bins, 
#                                            range=[-v_range,v_range],
#                                            normed=1, 
#                                            facecolor=colors[state],
#                                            weights = velo_data[:,1],
#                                            alpha=0.5, 
#                                            label='{}-bin{}'.format(state, cbin_id))
            plt.bar(centers, 
                    hist, 
                    align ='center', 
                    width = width,
                    alpha = 0.5,
                    color = colors[state],
                    label = '{}-bin{}'.format(state, cbin_id))
        plt.legend()
        plt.savefig('ana_trace_flux.velo.bin{}.pdf'.format(cbin_id))  
#     plt.legend()
#     plt.savefig('ana_trace_flux.velocities.pdf')   
#     plt.show()

#Analytical calculation of probabilites for the two particle system
if args.k >0:
    def PMF(k, R):
        r_0 = 15
        return  k*(R-r_0)**4 - 2*k* (R-r_0)**2 - 2*constants.kT *math.log(R)
    def pPMF(k,R):
        return math.exp( - PMF(k,R)/constants.kT)
    
    def pState(k, a,b):
        return scipy.integrate.quad(lambda r:pPMF(k,r), a, b)[0]
    
    probability_state_A_analytical = pState(args.k, state_A[0], state_A[1]) / pState(args.k, 1e-3, 1e3 )
    probability_state_B_analytical = pState(args.k, state_B[0], state_B[1]) / pState(args.k, 1e-3, 1e3 )

# Output
block_size = args.bs

def rate(datapoints):
    """
    calculates rate from mean(fluxes)/mean(probabilities)
    """
    return np.mean(datapoints[:,0]) / np.mean(datapoints[:,1])
############################
#      Pretty print        #
############################
sys.stderr.write('- calculating rates with errors\n')
sys.stdout.write("States:  A=[{},{}]  B=[{},{}]\n".format(state_A[0], state_A[1], state_B[0], state_B[1]))
# State A 
sys.stdout.write(" State A -> B:\n")
sys.stdout.flush()
if args.rates_only == False:
    block_bootstrap_flux_into_B =  f.block_bootstrap(flux_into_B[b:e], 
                                                     np.mean, 
                                                     block_size, 
                                                     number_of_samples = args.bs_samples)
    sys.stdout.write("    Flux(A->B):  {:5.4e}  CI: [{:5.4e}, {:5.4e}]\n".format(block_bootstrap_flux_into_B[0], 
                                                                   block_bootstrap_flux_into_B[1][0],
                                                                   block_bootstrap_flux_into_B[1][1]))
    sys.stdout.flush()
    block_bootstrap_prop_A = f.block_bootstrap(probability_state_A[b:e], 
                                               np.mean, 
                                               block_size, 
                                               number_of_samples = args.bs_samples)
    sys.stdout.write("    P(A):        {:5.4e}  CI: [{:5.4e}, {:5.4e}]\n".format(block_bootstrap_prop_A[0],
                                                                   block_bootstrap_prop_A[1][0],
                                                                   block_bootstrap_prop_A[1][1]))
    sys.stdout.flush()
    block_bootstrap_prop_from_A = f.block_bootstrap(probability_from_A[b:e], 
                                                    np.mean, 
                                                    block_size, 
                                                    number_of_samples = args.bs_samples)
    sys.stdout.write("    P(from A):   {:5.4e}  CI: [{:5.4e}, {:5.4e}]\n".format(block_bootstrap_prop_from_A[0],
                                                                   block_bootstrap_prop_from_A[1][0],
                                                                   block_bootstrap_prop_from_A[1][1]))
    sys.stdout.flush()

    if args.k >0:
        sys.stdout.write("    P(A) analytical: {:5.4e}\n".format(probability_state_A_analytical))
	sys.stdout.flush()

flux_prob_pair = []
for i in range(b,e):
    flux_prob_pair.append([flux_into_B[i], probability_state_A[i]])
flux_prob_pair = np.array(flux_prob_pair)

#print flux_prob_pair[0:10], np.mean(flux_prob_pair[0:10,0])
#print flux_into_B[b:b+10], probability_state_A[b:b+10], np.mean(flux_into_B[b:b+10]), np.mean(probability_state_A[b:b+10])
#print np.array(flux_prob_pair)
block_bootstrap_rate_into_B = f.block_bootstrap(flux_prob_pair, 
                                                rate, 
                                                block_size, 
                                                number_of_samples = args.bs_samples)
sys.stdout.write("    k {:5.4e} {:5.4e} {:5.4e}\n".format(block_bootstrap_rate_into_B[0], 
                                             block_bootstrap_rate_into_B[1][0],
                                             block_bootstrap_rate_into_B[1][1]))    
sys.stdout.flush()

# MFPT rate
flux_histprob_pair = []
for i in range(b,e):
    flux_histprob_pair.append([flux_into_B[i], probability_from_A[i]])
flux_histprob_pair = np.array(flux_histprob_pair)

block_bootstrap_mfpt_into_B = f.block_bootstrap(flux_histprob_pair, 
                                                rate, 
                                                block_size, 
                                                number_of_samples = args.bs_samples)
print("1/(MFPT) {:5.4e} {:5.4e} {:5.4e}".format(block_bootstrap_mfpt_into_B[0], 
                                             block_bootstrap_mfpt_into_B[1][0],
                                             block_bootstrap_mfpt_into_B[1][1]))


# State B
sys.stdout.write(" State B -> A:\n")
if args.rates_only == False:
    block_bootstrap_flux_into_A =  f.block_bootstrap(flux_into_A[b:e], 
                                                     np.mean, 
                                                     block_size, 
                                                     number_of_samples = args.bs_samples)
    sys.stdout.write("    Flux(B->A):  {:5.4e}  CI: [{:5.4e}, {:5.4e}]\n".format(block_bootstrap_flux_into_A[0], 
                                                                   block_bootstrap_flux_into_A[1][0],
                                                                   block_bootstrap_flux_into_A[1][1]))
    sys.stdout.flush()
    block_bootstrap_prop_B = f.block_bootstrap(probability_state_B[b:e], 
                                               np.mean, 
                                               block_size, 
                                               number_of_samples = args.bs_samples)
    sys.stdout.write("    P(B):        {:5.4e}  CI: [{:5.4e}, {:5.4e}]\n".format(block_bootstrap_prop_B[0],
                                                                   block_bootstrap_prop_B[1][0],
                                                                   block_bootstrap_prop_B[1][1]))
    sys.stdout.flush()
    block_bootstrap_prop_from_B = f.block_bootstrap(probability_from_B[b:e], 
                                                    np.mean, 
                                                    block_size, 
                                                    number_of_samples = args.bs_samples)
    sys.stdout.write("    P(from B):   {:5.4e}  CI: [{:5.4e}, {:5.4e}]\n".format(block_bootstrap_prop_from_B[0],
                                                                   block_bootstrap_prop_from_B[1][0],
                                                                   block_bootstrap_prop_from_B[1][1]))
    sys.stdout.flush()
    if args.k >0:
        sys.stdout.write("    P(B) analytical: {:5.4e}\n".format(probability_state_B_analytical))
        
flux_prob_pair = []
for i in range(b,e):
    flux_prob_pair.append([flux_into_A[i], probability_state_B[i]])
flux_prob_pair = np.array(flux_prob_pair)

block_bootstrap_rate_into_A = f.block_bootstrap(flux_prob_pair, 
                                                rate, 
                                                block_size,
                                                number_of_samples = args.bs_samples)
sys.stdout.write("    k {:5.4e} {:5.4e} {:5.4e}\n".format(block_bootstrap_rate_into_A[0], 
                                                           block_bootstrap_rate_into_A[1][0],
                                                           block_bootstrap_rate_into_A[1][1]))

# MFPT rate
flux_histprob_pair = []
for i in range(b,e):
    flux_histprob_pair.append([flux_into_A[i], probability_from_B[i]])
flux_histprob_pair = np.array(flux_histprob_pair)

block_bootstrap_mfpt_into_A = f.block_bootstrap(flux_histprob_pair, 
                                                rate, 
                                                block_size, 
                                                number_of_samples = args.bs_samples)
print("1/(MFPT) {:5.4e} {:5.4e} {:5.4e}".format(block_bootstrap_mfpt_into_A[0], 
                                             block_bootstrap_mfpt_into_A[1][0],
                                             block_bootstrap_mfpt_into_A[1][1]))  

########################################

