#!/usr/bin/python2
from __future__ import print_function
from logger import Logger
import analysis_operations
import argparse 
import numpy
import math
import matplotlib.pyplot as plt

parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")
parser.add_argument('-i', '--bin_index', dest="bin_index",
                    type=int, default=0,
                    help="Bin index.") 
parser.add_argument('-r', '--reweighting-range', dest="reweighting_range",
                    type=float, default=0.0,
                    help="reweighting range for mean value plotting.") 
parser.add_argument('-o', '--output', dest="output_filename",
                    type=str, default='OutrateEvol',
                    help="Output Filename.")
parser.add_argument('--ln', dest="ln", action="store_true",
                    default=False, help="print ln of rates")
                    
# Initialize
print('\033[1mAnalyzing Rates.\033[0m')      
args = parser.parse_args()


#get the actual Iteration from logger module
logger = Logger(args.logfile, APPEND = True)
iterations = logger.loadIterations(bCheckFiles=False)
logger.close()

#functions:
outrates = numpy.zeros([len(iterations), iterations[-1].getNumberOfBins()], float)
for i in range(0,len(outrates[:,0])):
    for j in range(0,len(outrates[0,:])):
        outrates[i,j] = 'inf'
        
for i in range(0, len(iterations)):
    if iterations[i].getNumberOfBins() >= args.bin_index+1:
        if args.reweighting_range > 0.0:
            iteration_range = int((iterations[i].getId() + 1) * args.reweighting_range)
        else:
            iteration_range = 0
        rates_tmp = analysis_operations.meanRateMatrix(iterations,i-iteration_range,i)[args.bin_index,:]
        for j in range(0, len(rates_tmp)):
            if rates_tmp[j] > 0.0:
                if args.ln == True:
                    outrates[i,j]  = math.log(rates_tmp[j])
                else:
                    outrates[i,j]  = rates_tmp[j]   
        
#Save to file
#header_line = 'Free Energy at: Bin, Iteration'            
#numpy.savetxt(md_module.workdir+args.output_path, bin_probabilities[:,:,1], header = header_line)
    
#Plot
fig=plt.figure(figsize=(5,5))
plt.title('rates out of bin #' + str(args.bin_index))
plt.xlabel('into bin #')
plt.ylabel('iteration')
plt.imshow(outrates, interpolation='none',origin='lower')
cbar=plt.colorbar()
plt.jet()
if args.ln == True:
    cbar.set_label('ln(rate per dt)')
else:
    cbar.set_label('rate per dt')
    
plt.savefig(args.output_filename+'.pdf',format='pdf',dpi=300)   

#print('\n Output written to: ' + args.output_path)
