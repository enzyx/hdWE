#!/usr/bin/python2
"""
Calculates the PMF along an arbitrary coordinate from 
the data of a hdWE run. Last iteration is including.
"""
import numpy
import sys
from logger import Logger
from math import log
import matplotlib.pyplot as plt
import constants

# Compatibility mode for python2.6
has_argparse = False
try:
    import argparse  
    has_argparse = True  
except ImportError:
    import optparse  #Python 2.6

###### Parse command line ###### 
if has_argparse:
    parser =argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
else:
    parser = optparse.OptionParser()
    parser.add_argument = parser.add_option

parser.add_argument('-d', '--dir', type=str, 
                    dest="work_dir", required=True, metavar="DIR",
                    help="The working direcory")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    required=False, type=int, default=0,
                    help="First iteration to use for PMF calculation.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    required=False, type=int, default=0,
                    help="Last iteration to to use for PMF calculation.")  
parser.add_argument('-o', '--output', dest="output_path", 
                    required=False, type=str, default='BinProbabilityEvolution',
                    help="Output filename")  
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    required=True, default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")#
parser.add_argument("--probability", dest="probability", action="store_true",
                    default=False)

print('\033[1mCalculating Bin Free Energies\033[0m (Free Energy is given in kcal/mol at 298K).')   
                  
args = parser.parse_args()

#get the actual Iteration from logger module
logger = Logger(args.work_dir+args.logfile)
iterations = logger.load_iterations()
logger.close()

#iteration I = iterations[I - 1]
args.first_iteration -= 1
args.last_iteration  -= 1

#set default last_iteration value
if args.last_iteration < 1:
    args.last_iteration = len(iterations) - 1


#initialize bin probability evolution array with size of last frame number_of_bins
n_iterations = args.last_iteration - args.first_iteration + 1
bin_probabilities = numpy.zeros([n_iterations, iterations[args.last_iteration].getNumberOfBins(),2], float)
for i in range(0,len(bin_probabilities[:,0,0])):
    for j in range(0,len(bin_probabilities[0,:,0])):
        for k in range(0,2):
            bin_probabilities[i,j,k] = 'Inf'

for i in range(args.first_iteration,args.last_iteration + 1):
    sys.stdout.write(' Processing iteration ' + str(iterations[i].getId()).zfill(5) +  \
                     ' / ' + str(args.first_iteration+1).zfill(5) + '-' + str(args.last_iteration+1).zfill(5) + '\r')
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
    numpy.savetxt(args.work_dir+args.output_path, bin_probabilities[:,:,1], header = header_line)
    
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
plt.savefig(args.work_dir+args.output_path+'.png',format='png',dpi=300)   

print('\n Output written to: ' + args.output_path)  


            


