#!/usr/bin/python2
"""
Calculates the PMF along an arbitrary coordinate from 
the data of a hdWE run. Last iteration is including.
"""
from __future__ import print_function
import numpy
import sys
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
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration)

#initialize bin probability evolution array with size of last frame number_of_bins
n_iterations = len(iterations)
bin_probabilities = numpy.zeros([n_iterations, iterations[-1].getNumberOfBins(),2], float)

for i in range(0,len(bin_probabilities[:,0,0])):
    for j in range(0,len(bin_probabilities[0,:,0])):
        for k in range(0,2):
            bin_probabilities[i,j,k] = 'Inf'



for i in range(n_iterations):
    for j in range(iterations[i].getNumberOfBins()):
        bin_probabilities[i,j,0] = iterations[i].bins[j].getProbability()
        if bin_probabilities[i,j,0] > 0.0:
            bin_probabilities[i,j,1] = - constants.kT * log(bin_probabilities[i,j,0]) 
        else:
            bin_probabilities[i,j,1] = 'Inf'
    minimum_free_energy = min(bin_probabilities[i,:,1])
    for j in range(0,iterations[i].getNumberOfBins()):
        bin_probabilities[i,j,1] -= minimum_free_energy        

# resort bins to represent the more meaningful coordinate sorting
resort_indices = binIdToCoordinateId(iterations[-1])
tmp_bin_probabilities = numpy.zeros([n_iterations, iterations[-1].getNumberOfBins(),2], float)
for i in range(bin_probabilities.shape[0]):
    for j in range(bin_probabilities.shape[1]):
        tmp_bin_probabilities[i,resort_indices[j],0] = bin_probabilities[i,j,0]
        tmp_bin_probabilities[i,resort_indices[j],1] = bin_probabilities[i,j,1]
bin_probabilities = tmp_bin_probabilities 

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


            


