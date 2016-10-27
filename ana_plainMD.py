#!/usr/bin/python2.7
#
# This file is part of hdWE. 
# Copyright (C) 2016 Manuel Luitz <manuel.luitz@tum.de>
# Copyright (C) 2016 Rainer Bomblies <r.bomblies@tum.de>
# Copyright (C) 2016 Fabian Zeller
#
# hdWE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# hdWE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with hdWE. If not, see <http://www.gnu.org/licenses/>.
# 
from __future__ import print_function
import sys
import argparse
import numpy
import lib.constants as constants
import scipy.stats
import lib.functions_ana_plainMD as functions_ana_plainMD
from math import log
from lib import functions_ana_general
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

sys.stderr.write('\033[1mAnalyzing plain MD trajectory\033[0m\n') 

args = parser.parse_args()

# Load coordinates
sys.stderr.write(' Loading datafile.\n')
coordinates_tmp = numpy.loadtxt(args.cpptraj_output, usecols=(args.use_column,) )
if args.end_frame == -1:
    args.end_frame = len(coordinates_tmp) - 1
coordinates = coordinates_tmp[args.begin_frame:args.end_frame]

sys.stderr.write(' Loaded frames: {:d}\n'.format(len(coordinates)))

######## PMF Calculation ######################################################
sys.stderr.write(' Calculating PMF.\n')
#Calculate the weighted histogram and PMF
pmf = functions_ana_plainMD.PMF(coordinates, args.number_of_bins)
pmf_header_line = 'Coordinate Value, Free energy, Probability'
numpy.savetxt(args.pmf_output, pmf, header = pmf_header_line)

######## Rate Calculation #####################################################
if not(args.state_A==None and args.state_B==None):
    sys.stderr.write(' Calculating Rates between States.\n')
    
    state_A = args.state_A
    state_B = args.state_B
    state_A = numpy.sort(state_A)
    state_B = numpy.sort(state_B)

    transitions_into_B =  functions_ana_plainMD.transitions_from_coordinates(coordinates_tmp, state_A, state_B)
    transitions_into_A =  functions_ana_plainMD.transitions_from_coordinates(coordinates_tmp, state_B, state_A)
    fpt_into_B         = functions_ana_plainMD.mfpt_from_coordinates(coordinates_tmp, state_A, state_B)
    fpt_into_A         = functions_ana_plainMD.mfpt_from_coordinates(coordinates_tmp, state_B, state_A)
    N_transitions_into_B = len(transitions_into_B)    
    N_transitions_into_A = len(transitions_into_A)         
   
    # Cumulative rates
    rates_into_B_cum = 1./numpy.array(functions_ana_general.cumulative_mean(transitions_into_B))    
    rates_into_A_cum = 1./numpy.array(functions_ana_general.cumulative_mean(transitions_into_A))
    fpt_into_B_cum   = 1./numpy.array(functions_ana_general.cumulative_mean(fpt_into_B))    
    fpt_into_A_cum   = 1./numpy.array(functions_ana_general.cumulative_mean(fpt_into_A))
      
    numpy.savetxt('ana_plainMD.B.cum', rates_into_B_cum)
    numpy.savetxt('ana_plainMD.A.cum', rates_into_A_cum)
    numpy.savetxt('ana_plainMD.mfpt_into_B.cum', fpt_into_B_cum)
    numpy.savetxt('ana_plainMD.mfpt_into_A.cum', fpt_into_A_cum)
    
                   
    # error with our bootstrapping
    def gustav_hilfsfunktion(inlist):
        return 1./numpy.mean(inlist)
    rates_into_A = functions_ana_general.block_bootstrap(transitions_into_A, gustav_hilfsfunktion, 1, number_of_samples=10000, alpha=0.05)
    rates_into_B = functions_ana_general.block_bootstrap(transitions_into_B, gustav_hilfsfunktion, 1, number_of_samples=10000, alpha=0.05)
    mfpt_into_A  = functions_ana_general.block_bootstrap(fpt_into_A, gustav_hilfsfunktion, 1, number_of_samples=10000, alpha=0.05)
    mfpt_into_B  = functions_ana_general.block_bootstrap(fpt_into_B, gustav_hilfsfunktion, 1, number_of_samples=10000, alpha=0.05)



    ####### Output
    print('')
    print(' State A ---> State B')
    print('    Number of transitions:      {:9d}'.format(N_transitions_into_B))
    print('')
    print('    Mean Rate:                  {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          #mean = meanCI_rates_into_B[0],
          mean = rates_into_A[0],
          lCI  = rates_into_A[1][0],
          uCI  = rates_into_A[1][1]))
    print('    MFPT Rate:                  {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          #mean = meanCI_rates_into_B[0],
          mean = mfpt_into_A[0],
          lCI  = mfpt_into_A[1][0],
          uCI  = mfpt_into_A[1][1]))

    print('')       
    print(' State B ---> State A')
    print('    Number of transitions:      {:9d}'.format(N_transitions_into_A))
    print('')
    print('    Mean Rate:                  {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          mean = rates_into_B[0],
          lCI  = rates_into_B[1][0],
          uCI  = rates_into_B[1][1]))
    print('    MFPT Rate:                  {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
          #mean = meanCI_rates_into_B[0],
          mean = mfpt_into_B[0],
          lCI  = mfpt_into_B[1][0],
          uCI  = mfpt_into_B[1][1]))

    print('')
        

