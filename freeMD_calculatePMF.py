import argparse
import numpy
import iteration
from math import log
from amber_module import MD_module

# Parse command line
parser = argparse.ArgumentParser(description=
    'Calculates the PMF based on a cpptraj output file generated from a free MD simulation.')
parser.add_argument('-f', '--first_frame', dest="first_frame",
                    required=False, type=int, default=0,
                    help="First frame to use for PMF calculation.")                    
parser.add_argument('-l', '--last_frame', dest="last_frame",
                    required=False, type=int, default=0,
                    help="Last frame to to use for PMF calculation.")  
parser.add_argument('-i', '--input', dest="input_path", 
                    required=True, type=str,
                    help="Output filename") 
parser.add_argument('-o', '--output', dest="output_path", 
                    required=False, type=str, default='freeMD_calculatePMF.output',
                    help="Output filename")  
parser.add_argument('-b', '--number_of_bins', dest="number_of_bins",
                    required=False, type=int, default=100, 
                    help="Number of bins used to calculate the probability histogram.")  

# Initialize
args = parser.parse_args()
kT = 0.598 # at 298K in kcal/mol 

# Load coordinates
coordinates_tmp = numpy.loadtxt(args.input_path, usecols=(1,) )
print len(coordinates_tmp)
if args.last_frame==0:
    args.last_frame = len(coordinates_tmp) - 1

coordinates = numpy.zeros([args.last_frame - args.first_frame])
coordinates = coordinates_tmp[args.first_frame:args.last_frame]
print coordinates

#Calculate the weighted histogram and PMF     
#Setup variables
hist_min =  min(coordinates)
hist_max =  max(coordinates)
print hist_max
dcoord   =  1.0 * (hist_max - hist_min ) / args.number_of_bins
hist     =  numpy.zeros([args.number_of_bins,2], float)
#Sort coords into histogramm
for i in range(0,len(coordinates)):
    index       = int( (coordinates[i] - hist_min) / dcoord )
    #maximum coord entry shall not be in an extra bin:
    if index==args.number_of_bins:
        index = index - 1
    hist[index,1] = hist[index,1] + 1
#Assign the bin positions and calculate free energy:
for i in range(0,args.number_of_bins):
    hist[i,0] = hist_min + i * dcoord
    if hist[i,1]>0:
        hist[i,1]  = - kT * log(hist[i,1])
    else:
        hist[i,1]  = 0

#Shift minimum to zero        
pmf_min = min(hist[:,1])
for i in range(0,args.number_of_bins):
    hist[i,1] = hist[i,1] - pmf_min

#Save PMF to file
header_line = 'Coordinate Value, Free energy'
numpy.savetxt(args.output_path, hist, header = header_line)


            


