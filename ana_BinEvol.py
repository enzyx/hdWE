#!/usr/bin/python2
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
"""
Calculates the PMF along an arbitrary coordinate from 
the data of a hdWE run. Last iteration is including.
"""
from __future__ import print_function
import numpy
from lib.logger import Logger
import lib.constants as constants
from math import log
import matplotlib.pyplot as plt
import argparse 
from lib.functions_ana_general import binIdToCoordinateId 

###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    default="hdWE-log", metavar="DIR",
                    help="The logfile for reading and writing")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    required=False, type=int, default=0,
                    help="First iteration to use for PMF calculation.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    required=False, type=int, default=-1,
                    help="Last iteration to to use for PMF calculation.")  
parser.add_argument('-o', '--output', dest="output_path", 
                    required=True, type=str, default='BinProbabilityEvolution',
                    help="Output filename")  
parser.add_argument("--probability", dest="probability", action="store_true",
                    default=False)
                  
# Initialize
print('\033[1mCalculating Bin Free Energies\033[0m (Free Energy is given in kcal/mol at 298K).')
args = parser.parse_args()
last_iteration = args.last_iteration 
first_iteration = args.first_iteration
logger = Logger(args.logdir)
if last_iteration < 0:
    last_iteration = logger.getLastIterationId()
    
n_iterations = last_iteration - first_iteration + 1 
    
#initialize bin probability evolution array with size of last frame number_of_bins
iteration = logger.loadIterations(logger.getLastIterationId(), logger.getLastIterationId())[0]
n_bins = iteration.getNumberOfBins()
bin_probabilities = numpy.zeros([n_iterations, iteration.getNumberOfBins(),2], float)
for i in range(0,len(bin_probabilities[:,0,0])):
    for j in range(0,len(bin_probabilities[0,:,0])):
        for k in range(0,2):
            bin_probabilities[i,j,k] = 0





for i in range(first_iteration, last_iteration):
    print(i)
    iteration = logger.loadIterations(i, i)[0]
    i = i - first_iteration
    for j in range(iteration.getNumberOfBins()):
        bin_probabilities[i,j,0] = iteration.bins[j].getProbability()
        if bin_probabilities[i,j,0] > 0.0:
            bin_probabilities[i,j,1] = - constants.kT * log(bin_probabilities[i,j,0]) 
        else:
            bin_probabilities[i,j,1] = 'Inf'
    minimum_free_energy = min(bin_probabilities[i,:,1])
    for j in range(0,iteration.getNumberOfBins()):
        bin_probabilities[i,j,1] -= minimum_free_energy        

# resort bins to represent the more meaningful coordinate sorting
resort_indices = binIdToCoordinateId(iteration)
print(resort_indices)
tmp_bin_probabilities = numpy.zeros([n_iterations, n_bins,2], float)
for i in range(bin_probabilities.shape[0]):
    for j in range(bin_probabilities.shape[1]):
        tmp_bin_probabilities[i,resort_indices[j],0] = bin_probabilities[i,j,0]
        tmp_bin_probabilities[i,resort_indices[j],1] = bin_probabilities[i,j,1]
bin_probabilities = tmp_bin_probabilities 
mean_bin_free_energies = numpy.zeros([n_bins,2], float)
for j in range(0,n_bins):
    if j == 0:
        mean_bin_free_energies[j,0] = iteration.boundaries[0][0]
    elif j == (n_bins - 1):
        mean_bin_free_energies[j,0] = iteration.boundaries[0][j-1]
    else:
        mean_bin_free_energies[j,0] = (iteration.boundaries[0][j-1] + iteration.boundaries[0][j] ) / 2
    mean_bin_free_energies[j,1] = numpy.sum(bin_probabilities[:,j,0])

#for j in range(0,len(mean_bin_free_energies)):
#        if mean_bin_free_energies[j,1] > 0.0:
#            mean_bin_free_energies[j,1] = - constants.kT * log(mean_bin_free_energies[j,1]) 
#
#minimum_free_energy = min(mean_bin_free_energies[:,1])
#
#for j in range(0,len(mean_bin_free_energies)):
#    mean_bin_free_energies[j,1] -= minimum_free_energy 
    
    
numpy.savetxt(args.output_path+'mean', mean_bin_free_energies)   
#Save to file
if args.probability == True:
    header_line = 'Probability at: Bin, Iteration'            
    numpy.savetxt(args.output_path, bin_probabilities[:,:,0], header = header_line)
else:
    header_line = 'Free Energy at: Bin, Iteration'            
    numpy.savetxt(args.output_path, bin_probabilities[:,:,1], header = header_line)
 
# Plot as png
fig=plt.figure(figsize=(5,4))
plt.xlabel('# iteration')
plt.ylabel('# bin')
if args.probability == True:
    plt.imshow(numpy.transpose(bin_probabilities[:,:,0]), interpolation='none',origin='lower')
    cbar=plt.colorbar()
    plt.jet()
    cbar.set_label('Probability')
else:
    plt.imshow(numpy.transpose(bin_probabilities[:,:,1]), interpolation='none',origin='lower')
    #cbar = plt.colorbar(fraction=0.1)
    plt.jet()
    #cbar.set_label('Free Energy in kT')
plt.savefig(args.output_path+'.png',format='png',dpi=300)   

print('\n Output written to: ' + args.output_path)


            


