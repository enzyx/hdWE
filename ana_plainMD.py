#!/usr/bin/python2.7
from __future__ import print_function
import argparse
import numpy
import lib.constants as constants
import scipy.stats
import lib.functions_ana_plainMD as functions_ana_plainMD
from math import log
#import scikits.bootstrap

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
parser.add_argument('--output-pmf', dest="pmf_output", 
                    required=False, type=str, default='plainMD.pmf',
                    help="PMF output filename")  
parser.add_argument('--output-transition-time-distribution-into-A', dest="distr_into_A_output", 
                    required=False, type=str, default='plainMD.trans_time_distr_into_A',
                    help="Transition time into state A distribution output filename") 
parser.add_argument('--output-transition-time-distribution-into-B', dest="distr_into_B_output", 
                    required=False, type=str, default='plainMD.trans_time_distr_into_B',
                    help="Transition time into state B distribution output filename")  
parser.add_argument('-N', '--number_of_bins', dest="number_of_bins",
                    required=False, type=int, default=200, 
                    help="Number of bins used to calculate the probability histogram.")  
parser.add_argument('--state-A', dest="state_A",
                    required=False, type=float, nargs=2, 
                    help="Boundaries of the start state for rate calculation.")  
parser.add_argument('--state-B', dest="state_B",
                    required=False, type=float, nargs=2, 
                    help="Boundaries of the end state for rate calculation.") 
parser.add_argument('--N-tau-per-WE-iteration', dest="N_tau_per_it",
                    required=False, type=int, default = 0,
                    help="Number of calculated taus per WE iteration for cumulative plot comparable to WE run.") 
parser.add_argument('-c', '--input-column', dest="use_column",
                    required=False, type=int, default=1,
                    help="The column with coordinate data in input file.")
parser.add_argument('--passage-times', nargs='?', default=False, const='fpt.dat', 
                    dest='write_mfpt', metavar='FILE',
                    help="Drop first passage times to (this) file.")

######## Initialize ###########################################################

print('\033[1mAnalyzing plain MD trajectory\033[0m') 

args = parser.parse_args()

# Load coordinates
coordinates_tmp = numpy.loadtxt(args.cpptraj_output, usecols=(args.use_column,) )
if args.end_frame == -1:
    args.end_frame = len(coordinates_tmp) - 1
coordinates = coordinates_tmp[args.begin_frame:args.end_frame]

print(' Loaded frames: {:d}'.format(len(coordinates)))

######## PMF Calculation ######################################################
print(' Calculating PMF.')
#Calculate the weighted histogram and PMF
pmf = functions_ana_plainMD.PMF(coordinates, args.number_of_bins)
pmf_header_line = 'Coordinate Value, Free energy, Probability'
numpy.savetxt(args.pmf_output, pmf, header = pmf_header_line)

######## Rate Calculation #####################################################
if not(args.state_A==None and args.state_B==None):
    print(' Calculating Rates between States.')
    
    state_A = args.state_A
    state_B = args.state_B
    state_A = numpy.sort(state_A)
    state_B = numpy.sort(state_B)

    ####### Cumulative plot
    if args.N_tau_per_it > 0:
        fout = open('cumulative_flux.dat', 'w')
        fout.write("# F->A           P_A       F->B            P_B           \n")
        fout.close() 
        for i in range(1, int(len(coordinates) / args.N_tau_per_it) ):  
            print(i)  
            coordinates_tmp = coordinates[:i*args.N_tau_per_it]
            first_passage_times_into_B, transition_times_into_B =  functions_ana_plainMD.transitions_from_coordinates(coordinates_tmp, state_A, state_B)
            first_passage_times_into_A, transition_times_into_A =  functions_ana_plainMD.transitions_from_coordinates(coordinates_tmp, state_B, state_A)
            N_transitions_into_B = len(transition_times_into_B)    
            N_transitions_into_A = len(transition_times_into_A)
            residence_times_in_A = functions_ana_plainMD.residence_times_from_coordinates(coordinates_tmp, state_A)    
            residence_times_in_B = functions_ana_plainMD.residence_times_from_coordinates(coordinates_tmp, state_B) 
            flux_into_B          = 1.0*N_transitions_into_B / len(coordinates_tmp)
            flux_into_A          = 1.0*N_transitions_into_A / len(coordinates_tmp)
            probability_in_A     = 1.0*numpy.sum(residence_times_in_A) / len(coordinates_tmp)
            probability_in_B     = 1.0*numpy.sum(residence_times_in_B) / len(coordinates_tmp)
            # Data time series
            fout = open('cumulative_flux.dat', 'a')
            fout.write("{:8.7e}   {:8.7e}    {:8.7e}   {:8.7e}\n".format(
                                 flux_into_A, 
                                 probability_in_A,  
                                 flux_into_B, 
                                 probability_in_B) )
            fout.close()    

    ####### times
    first_passage_times_into_B, transition_times_into_B =  functions_ana_plainMD.transitions_from_coordinates(coordinates, state_A, state_B)
    first_passage_times_into_A, transition_times_into_A =  functions_ana_plainMD.transitions_from_coordinates(coordinates, state_B, state_A)
    
    N_transitions_into_B = len(transition_times_into_B)    
    N_transitions_into_A = len(transition_times_into_A)
    
    residence_times_in_A = functions_ana_plainMD.residence_times_from_coordinates(coordinates, state_A)    
    residence_times_in_B = functions_ana_plainMD.residence_times_from_coordinates(coordinates, state_B)  
   
   
    # write first passage times to file
    if args.write_mfpt:
        fpt_into_A = numpy.array(first_passage_times_into_A)
        fpt_into_B = numpy.array(first_passage_times_into_B)
        # get them to equal size
        while len(fpt_into_A) < len(fpt_into_B):
            fpt_into_B = fpt_into_B[:-1]
        while len(fpt_into_A) > len(fpt_into_B):
            fpt_into_A = fpt_into_A[:-1]            
        mfpt = numpy.transpose(numpy.matrix([fpt_into_A, fpt_into_B]))
        #print (mfpt)
        numpy.savetxt(args.write_mfpt, mfpt)
   
    ####### Transition time distribution
    distr_header_line = 'Transition time in tau, Probability'
#     if len(transition_times_into_B) > 2:
#         transition_into_B_distribution = functions_ana_general.histogram(transition_times_into_B, 100)
#         numpy.savetxt(args.distr_into_B_output, transition_into_B_distribution, header = distr_header_line)
#     else:
#         print('No transitions to state B.')        
#     
#     if len(transition_times_into_A) > 2:
#         transition_into_A_distribution = functions_ana_general.histogram(transition_times_into_A, 100)
#         numpy.savetxt(args.distr_into_A_output, transition_into_A_distribution, header = distr_header_line)
#     else:
#         print('No transitions to state A.') 
#                  
    ####### 95% Confidence intervals
    # First Passage Times                 
    meanCI_first_passage_times_into_B = scipy.stats.bayes_mvs(first_passage_times_into_B, 0.95)[0]
    meanCI_first_passage_times_into_A = scipy.stats.bayes_mvs(first_passage_times_into_A, 0.95)[0]
    # Transition Times
    meanCI_transition_times_into_B    = scipy.stats.bayes_mvs(transition_times_into_B, 0.95)[0]
    meanCI_transition_times_into_A    = scipy.stats.bayes_mvs(transition_times_into_A, 0.95)[0]
    # Residence Times
    # Workaround because of bootstrap memory error for large coordinate sets:
    #sumCI_residence_times_in_A          = [numpy.sum(residence_times_in_A)]
    #sumCI_residence_times_in_A.append(scikits.bootstrap.ci(residence_times_in_A, numpy.sum, alpha = 0.05))
    sumCI_residence_times_in_A        = scipy.stats.bayes_mvs(residence_times_in_A, 0.95)[0]
    sumCI_residence_times_in_A        = (sumCI_residence_times_in_A[0] * len(residence_times_in_A), (sumCI_residence_times_in_A[1][0] * len(residence_times_in_A),sumCI_residence_times_in_A[1][1] * len(residence_times_in_A)))
    #sumCI_residence_times_in_B          = [numpy.sum(residence_times_in_B)]
    #sumCI_residence_times_in_B.append(scikits.bootstrap.ci(residence_times_in_B, numpy.sum, alpha = 0.05))
    sumCI_residence_times_in_B        = scipy.stats.bayes_mvs(residence_times_in_B, 0.95)[0]
    sumCI_residence_times_in_B        = (sumCI_residence_times_in_B[0] * len(residence_times_in_B), (sumCI_residence_times_in_B[1][0] * len(residence_times_in_B),sumCI_residence_times_in_B[1][1] * len(residence_times_in_B)))
    # Rates
    meanCI_rates_into_B               = ( functions_ana_plainMD.rate(N_transitions_into_B, sumCI_residence_times_in_A[0]), ( functions_ana_plainMD.rate(N_transitions_into_B, sumCI_residence_times_in_A[1][1]), functions_ana_plainMD.rate(N_transitions_into_B, sumCI_residence_times_in_A[1][0]) )) 
    meanCI_rates_into_A               = ( functions_ana_plainMD.rate(N_transitions_into_A, sumCI_residence_times_in_B[0]), ( functions_ana_plainMD.rate(N_transitions_into_A, sumCI_residence_times_in_B[1][1]), functions_ana_plainMD.rate(N_transitions_into_A, sumCI_residence_times_in_B[1][0]) )) 

    ####### Compare free energy difference between states
    # from PMF
    F_A = functions_ana_plainMD.F_of_state_from_hist(pmf, state_A)
    F_B = functions_ana_plainMD.F_of_state_from_hist(pmf, state_B)
    dF_pmf = F_B - F_A
    # from rates
    dF_rates = - constants.kT*log(meanCI_rates_into_B[0] / meanCI_rates_into_A[0])


    ####### Output
    print('')
    print(' State A ---> State B')
    print('    Number of transitions:      {:9d}'.format(N_transitions_into_B))
    print('    Mean First Passage time:    {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          mean = meanCI_first_passage_times_into_B[0],
          lCI  = meanCI_first_passage_times_into_B[1][0],
          uCI  = meanCI_first_passage_times_into_B[1][1]))
    print('    Mean Transition time:       {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          mean = meanCI_transition_times_into_B[0],
          lCI  = meanCI_transition_times_into_B[1][0],
          uCI  = meanCI_transition_times_into_B[1][1]))
    print('')
    print('    Total Residence time in A:  {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          mean = sumCI_residence_times_in_A[0],
          lCI  = sumCI_residence_times_in_A[1][0],
          uCI  = sumCI_residence_times_in_A[1][1]))
    print('    Mean Rate:                  {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          #mean = meanCI_rates_into_B[0],
          mean = meanCI_rates_into_B[0],
          lCI  = meanCI_rates_into_B[1][0],
          uCI  = meanCI_rates_into_B[1][1]))

    print('')       
    print(' State B ---> State A')
    print('    Number of transitions:      {:9d}'.format(N_transitions_into_A))
    print('    Mean First passage time:    {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          mean = meanCI_first_passage_times_into_A[0],
          lCI  = meanCI_first_passage_times_into_A[1][0],
          uCI  = meanCI_first_passage_times_into_A[1][1]))
    print('    Mean Transition time:       {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          mean = meanCI_transition_times_into_A[0],
          lCI  = meanCI_transition_times_into_A[1][0],
          uCI  = meanCI_transition_times_into_A[1][1]))
    print('')
    print('    Total Residence Time in A:  {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          mean = sumCI_residence_times_in_B[0],
          lCI  = sumCI_residence_times_in_B[1][0],
          uCI  = sumCI_residence_times_in_B[1][1]))
    print('    Mean Rate:                  {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          mean = meanCI_rates_into_A[0],
          lCI  = meanCI_rates_into_A[1][0],
          uCI  = meanCI_rates_into_A[1][1]))

    print('')       
    print('    dF from PMF:               {:5.3e}'.format(dF_pmf))
    print('    dF from Rates:             {:5.3e}'.format(dF_rates))
        

