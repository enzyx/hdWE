#!/usr/bin/python2
"""
Calculates the PMF along an arbitrary coordinate from 
the data of a hdWE run. Last iteration is including.
"""
from __future__ import print_function
import numpy
import sys
from logger import Logger
from math import log
import matplotlib.pyplot as plt
import constants
from amber_module import MD_module
import argparse 

###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-c', '--conf', dest="input_md_conf", 
                    type=str, required=True,
                    help="MD-Software configuration file")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    required=False, type=int, default=0,
                    help="First iteration to use for PMF calculation.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    required=False, type=int, default=-1,
                    help="Last iteration to to use for PMF calculation.")  
parser.add_argument('-o', '--output', dest="output_path", 
                    required=False, type=str, default='BinProbabilityEvolution',
                    help="Output filename")  
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    required=True, default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")#
parser.add_argument("--probability", dest="probability", action="store_true",
                    default=False)

                  
# Initialize
print('\033[1mCalculating Bin Free Energies\033[0m (Free Energy is given in kcal/mol at 298K).')   
args = parser.parse_args()
md_module = MD_module(args.input_md_conf, debug=False)

#get the actual Iteration from logger module
logger = Logger(args.logfile, append = True)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration)
logger.close()



#initialize bin probability evolution array with size of last frame number_of_bins
n_iterations = args.last_iteration - args.first_iteration + 1
bin_probabilities = numpy.zeros([n_iterations, iterations[args.last_iteration-1].getNumberOfBins(),2], float)
for i in range(0,len(bin_probabilities[:,0,0])):
    for j in range(0,len(bin_probabilities[0,:,0])):
        for k in range(0,2):
            bin_probabilities[i,j,k] = 'Inf'


for i in range(args.first_iteration,args.last_iteration):
    sys.stdout.write(' Processing iteration ' + str(i).zfill(5) +  \
                     ' / ' + str(args.first_iteration).zfill(5) + '-' + str(args.last_iteration).zfill(5) + '\r')
    for j in range(0,iterations[i].getNumberOfBins()):
        bin_probabilities[i,j,0] = iterations[i].bins[j].getProbability()
        if bin_probabilities[i,j,0] > 0.0:
            bin_probabilities[i,j,1] = - constants.kT * log(bin_probabilities[i,j,0]) 
        else:
            bin_probabilities[i,j,1] = 'Inf'
    minimum_free_energy = min(bin_probabilities[i,:,1])
    for j in range(0,iterations[i].getNumberOfBins()):
        bin_probabilities[i,j,1] -= minimum_free_energy        
      
#Save to file
if args.probability == True:
    header_line = 'Probability at: Bin, Iteration'            
    numpy.savetxt(args.work_dir+args.output_path, bin_probabilities[:,:,0], header = header_line)
else:
    header_line = 'Free Energy at: Bin, Iteration'            
    numpy.savetxt(md_module.workdir+args.output_path, bin_probabilities[:,:,1], header = header_line)
    
#Plot as png
fig=plt.figure(figsize=(5,5))
plt.xlabel('# bin')
plt.ylabel('# iteration')
if args.probability == True:
    plt.imshow(bin_probabilities[:,:,0], interpolation='none',origin='lower')
    cbar=plt.colorbar()
    plt.jet()
    cbar.set_label('Probability')
else:
    plt.imshow(bin_probabilities[:,:,1], interpolation='none',origin='lower')
    cbar=plt.colorbar()
    plt.jet()
    cbar.set_label('Free Energy in kT')
plt.savefig(md_module.workdir+args.output_path+'.png',format='png',dpi=300)   

print('\n Output written to: ' + args.output_path)  


            


