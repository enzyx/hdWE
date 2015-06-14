#!/usr/bin/python2.7
from __future__ import print_function
import argparse
import numpy
import constants
import scipy.stats
from math import log

# Parse command line
parser = argparse.ArgumentParser(description=
    'Calculates PMF and rates for a plain MD. The reaction coordinate cpptraj_output is required.')
parser.add_argument('-b', '--begin_frame', dest="begin_frame",
                    required=False, type=int, default=0,
                    help="First frame to use for PMF calculation.")                    
parser.add_argument('-e', '--end_frame', dest="end_frame",
                    required=False, type=int, default=-1,
                    help="Last frame to to use for PMF calculation.")  
parser.add_argument('-i', '--cpptraj_output', dest="cpptraj_output", 
                    type=str, 
                    help="File containing cpptraj output.")
parser.add_argument('-o', '--output', dest="output_path", 
                    required=False, type=str, default='freeMD_calculatePMF.output',
                    help="Output filename")  
parser.add_argument('-N', '--number_of_bins', dest="number_of_bins",
                    required=False, type=int, default=200, 
                    help="Number of bins used to calculate the probability histogram.")  
parser.add_argument('--state-A', dest="state_A",
                    required=False, type=float, nargs=2, 
                    help="Boundaries of the start state for rate calculation.")  
parser.add_argument('--state-B', dest="state_B",
                    required=False, type=float, nargs=2, 
                    help="Boundaries of the end state for rate calculation.")  


def transitions_from_coordinates(coordinates, start_state, end_state):
    transition_times    = []
    in_transition       = False
    transition_time_tmp = 0
    for coord in coordinates:
        if in_transition == False:
            # if in start state, start with transition time counting
            if (coord > start_state[0]) and (coord < start_state[1]):
                in_transition       = True
                transition_time_tmp = 0
        if in_transition == True:
            # if in transition, save transition time if arrived in end state,
            # otherwise add frame to transition time
            if (coord > end_state[0]) and (coord < end_state[1]):
                transition_times.append(transition_time_tmp)
                in_transition    = False
            else:
                transition_time_tmp += 1
    
    if len(transition_times) == 0:
        print(' No transitions found.')
    
    rates = numpy.zeros([len(transition_times)], float)
    for i in range(len(rates)):
        rates[i] = 1.0 / transition_times[i]        

    return transition_times, rates
    
def distribution_from_times(transition_times, N_bins):
    distr_min = min(transition_times)
    distr_max = max(transition_times)
    d         =  1.0 * (distr_max - distr_min ) / N_bins
    distr     =  numpy.zeros([N_bins,2], float)
    #Sort coords into histogramm
    for i in range(0,len(transition_times)):
        index       = int( (transition_times[i] - distr_min) / d )
        #maximum coord entry shall not be in an extra bin:
        if index==N_bins:
            index = index - 1
        distr[index,1] += 1
    #Assign the bin positions:
    for i in range(0,N_bins):
        distr[i,0] = distr_min + i * d
   
    return distr


# Initialize
args = parser.parse_args()

# Load coordinates
coordinates_tmp = numpy.loadtxt(args.cpptraj_output, usecols=(1,) )
if args.end_frame == -1:
    args.end_frame = len(coordinates_tmp) - 1

coordinates = numpy.zeros([args.end_frame - args.begin_frame])
coordinates = coordinates_tmp[args.begin_frame:args.end_frame]

######## PMF Calculation ######################################################

print('Calculating PMF.')


#Calculate the weighted histogram and PMF     
#Setup variables
hist_min =  min(coordinates)
hist_max =  max(coordinates)

dcoord   =  1.0 * (hist_max - hist_min ) / args.number_of_bins
hist     =  numpy.zeros([args.number_of_bins,3], float)
#Sort coords into histogramm
for i in range(0,len(coordinates)):
    index       = int( (coordinates[i] - hist_min) / dcoord )
    #maximum coord entry shall not be in an extra bin:
    if index==args.number_of_bins:
        index = index - 1
    hist[index,2] = hist[index,2] + 1
#Assign the bin positions and calculate free energy:
for i in range(0,args.number_of_bins):
    hist[i,0] = hist_min + i * dcoord
    if hist[i,2]>0:
        hist[i,1]  = - constants.kT * log(hist[i,2])
    else:
        hist[i,1]  = 'Inf'

#Shift minimum to zero        
pmf_min = min(hist[:,1])
for i in range(0,args.number_of_bins):
    hist[i,1] = hist[i,1] - pmf_min

#Save PMF to file
header_line = 'Coordinate Value, Free energy, Probability'
numpy.savetxt(args.output_path, hist, header = header_line)

######## Rate Calculation #####################################################

if not(args.state_A==None and args.state_B==None):
    print('Calculating rate from start state to end state.')
    
    state_A             = args.state_A
    state_B             = args.state_B

    transition_times_into_B, rates_into_B = transitions_from_coordinates(coordinates, state_A, state_B)
    transition_times_into_A, rates_into_A = transitions_from_coordinates(coordinates, state_B, state_A)
    
    transition_into_B_distribution = distribution_from_times(transition_times_into_B, 50)
    transition_into_A_distribution = distribution_from_times(transition_times_into_A, 50)
                    
    meanCI_transition_times_into_B = scipy.stats.bayes_mvs(transition_times_into_B, 0.95)[0]
    meanCI_transition_times_into_A = scipy.stats.bayes_mvs(transition_times_into_A, 0.95)[0]
    meanCI_rates_into_B            = scipy.stats.bayes_mvs(rates_into_B, 0.95)[0]
    meanCI_rates_into_A            = scipy.stats.bayes_mvs(rates_into_A, 0.95)[0]

    # print

    print(' State A ---> State B')
    print('    Number of transitions: {:9d}'.format(len(transition_times_into_B)))
    print('    Transition time:       {mean:5.3e}   CI: {lCI:5.3e} - {uCI:5.3e}'.format(
          mean = meanCI_transition_times_into_B[0],
          lCI  = meanCI_transition_times_into_B[1][0],
          uCI  = meanCI_transition_times_into_B[1][1]))
    print('    Rate:                  {mean:5.3e}   CI: {lCI:5.3e} - {uCI:5.3e}'.format(
          mean = meanCI_rates_into_B[0],
          lCI  = meanCI_rates_into_B[1][0],
          uCI  = meanCI_rates_into_B[1][1]))
          
    print(' State B ---> State A')
    print('    Number of transitions: {:9d}'.format(len(transition_times_into_A)))
    print('    Transition time:       {mean:5.3e}   CI: {lCI:5.3e} - {uCI:5.3e}'.format(
          mean = meanCI_transition_times_into_A[0],
          lCI  = meanCI_transition_times_into_A[1][0],
          uCI  = meanCI_transition_times_into_A[1][1]))
    print('    Rate:                  {mean:5.3e}   CI: {lCI:5.3e} - {uCI:5.3e}'.format(
          mean = meanCI_rates_into_A[0],
          lCI  = meanCI_rates_into_A[1][0],
          uCI  = meanCI_rates_into_A[1][1]))
          
    numpy.savetxt('distr_into_B', transition_into_B_distribution)
    numpy.savetxt('distr_into_A', transition_into_A_distribution)
