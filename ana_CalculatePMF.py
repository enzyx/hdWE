"""
Calculates the PMF along an arbitrary coordinate 
from the data of a hdWE run.
"""
from __future__ import print_function
import numpy
import constants
import sys
import segment
from logger import Logger
from math import log
from amber_module import MD_module
import argparse  


###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-c', '--conf', dest="input_md_conf", 
                    type=str, required=True,
                    help="MD-Software configuration file")
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0,
                    help="First iteration to use for PMF calculation.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Last iteration to to use for PMF calculation.")  
parser.add_argument('-o', '--output', dest="output_path", 
                    type=str, default='ana_calculatePMF.output',
                    help="Output filename")  
parser.add_argument('-N', '--number_of_bins', dest="number_of_bins",
                    type=int, default=100, 
                    help="Number of bins used to calculate the probability histogram.")  
parser.add_argument('-i', '--cpptraj_lines_file', dest="cpptraj_lines_file_path", 
                    type=str,   required=True,
                    help="File containig cpptraj syntax that defines the reaction coordinate.")
parser.add_argument('-y', dest="use_trajectory", action="store_true",
                    default=False, help="use conformations from the trajectories, not only the end conformation")                    
                    
# Initialize
print('\033[1mCalculating PMF\033[0m (Free Energy is given in kcal/mol at 298K).')      
args = parser.parse_args()
md_module = MD_module(args.input_md_conf, debug=False)

#get the actual Iteration from logger module
logger = Logger(args.logfile, APPEND = True)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration)
logger.close()

print(str(len(iterations)))

# Load cpptraj input file as one string with linebreaks and delete the last line break
try:
    cpptraj_lines_file=open(args.cpptraj_lines_file_path, 'r')
except:
    print('Error: could not open ' + args.cpptraj_lines_file_path)
cpptraj_lines=''
for line in cpptraj_lines_file:
    cpptraj_lines = cpptraj_lines + line
cpptraj_lines = cpptraj_lines[0:-1]
cpptraj_lines_file.close()

#Calculate the coordinate values and store them together with
#the trajectory probability into coordinates 
coordinates         = numpy.zeros([0,2])
coordinates_tmp     = numpy.zeros([1,2])
coordinates_tmp_tmp = numpy.zeros([1,2])

first_it = iterations[0].getId()
last_it  = iterations[-1].getId()
 
for iteration_loop in iterations:
    for bin_loop in iteration_loop:
        sys.stdout.write(' Processing iteration ' + str(iteration_loop.getId()).zfill(5) +  \
                         ' / ' + str(first_it).zfill(5) + '-' + str(last_it).zfill(5) + \
                         ', Bin ' +  str(bin_loop.getId()+1).zfill(5) + '/ ' + str(iteration_loop.getNumberOfBins()).zfill(5) + '\r' )
        sys.stdout.flush()  
        for segment_loop in bin_loop:
            coordinates_tmp_tmp = md_module.ana_calcCoordinateOfSegment(segment_loop, cpptraj_lines, args.use_trajectory)
            probability_tmp = segment_loop.getProbability()
            for i in range(0,len(coordinates_tmp_tmp)):
                coordinates_tmp[0,0] = coordinates_tmp_tmp[i]
                coordinates_tmp[0,1] = probability_tmp
                coordinates          = numpy.append(coordinates_tmp, coordinates, axis=0)

#Calculate the coordinate values of the bin reference structures
bin_coordinates     = numpy.zeros([0,1])
bin_coordinates_tmp = numpy.zeros([1,1])
 
for bin_loop in iterations[-1]:
    #create temporary segment to pass to md_module, because the bin reference segment
    #is not necessarilly within the loaded range of iterations
    segment_tmp = segment.Segment(probability = 0, parent_bin_id = 0, parent_segment_id = 0,
                  iteration_id = bin_loop.getReferenceIterationId(),
                  bin_id       = bin_loop.getReferenceBinId(),
                  segment_id   = bin_loop.getReferenceSegmentId() )
    bin_coordinates_tmp[0] = md_module.ana_calcCoordinateOfSegment(segment_tmp, cpptraj_lines, False)
    bin_coordinates        = numpy.append(bin_coordinates, bin_coordinates_tmp, axis=0) 

#Calculate the weighted histogram and PMF     
#Setup variables
hist_min =  min(coordinates[:,0])
hist_max =  max(coordinates[:,0])
dcoord   =  1.0 * (hist_max - hist_min ) / args.number_of_bins
hist     =  numpy.zeros([args.number_of_bins,5], float)
#Sort coords into histogramm
for i in range(0,len(coordinates[:,0])):
    index       = int( (coordinates[i,0] - hist_min) / dcoord )
    #maximum coord entry shall not be in an extra bin:
    if index==args.number_of_bins:
        index = index - 1
    hist[index,2] = hist[index,2] + coordinates[i,1]
    hist[index,3] = hist[index,3] + 1
#Sort bin coords into histogramm
for i in range(0,len(bin_coordinates[:,0])):
    index       = int( (bin_coordinates[i,0] - hist_min) / dcoord )
    #maximum bin coord entry shall not be in an extra bin:
    if index>=args.number_of_bins:
        index = args.number_of_bins - 1
    hist[index,4] = hist[index,4] + 1
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
header_line = 'Coordinate Value, Free Energy, Probability, Values per hist-Bin, hdWE Bins per hist-Bin'
numpy.savetxt(args.output_path, hist, header = header_line)

print('\n Output written to: ' + args.output_path)   

