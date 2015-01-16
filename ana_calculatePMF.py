import argparse
import numpy
import iteration
from math import log
from amber_module import MD_module

# Parse command line
parser = argparse.ArgumentParser(description=
    'Calculates the PMF along an arbitrary coordinate from the data of a hdWE run.')
parser.add_argument('-d', '--dir', type=str, 
                    dest="work_dir", required=True, metavar="DIR",
                    help="The working direcory")
parser.add_argument('-c', '--conf', dest="input_md_conf", 
                    required=True, type=str,
                    help="MD-Software configuration file")
parser.add_argument('-f', '--first_it', dest="first_iteration",
                    required=False, type=int, default=0,
                    help="First iteration to use for PMF calculation.")                    
parser.add_argument('-l', '--last_it', dest="last_iteration",
                    required=False, type=int, default=1,
                    help="Last iteration to to use for PMF calculation.")  
parser.add_argument('-o', '--output', dest="output_path", 
                    required=False, type=str, default='ana_calculatePMF.output',
                    help="Output filename")  
parser.add_argument('-b', '--number_of_bins', dest="number_of_bins",
                    required=False, type=int, default=100, 
                    help="Number of bins used to calculate the probability histogram.")  
parser.add_argument('-i', '--cpptraj_lines_file', dest="cpptraj_lines_file_path", 
                    required=True, type=str, 
                    help="File containig cpptraj syntax that defines the reaction coordinate.")
# Debug default False
# Free Energy unit default kcal/mol               
               
# Initialize
args = parser.parse_args()
md_module = MD_module(args.work_dir, args.input_md_conf, debug=False)
kT = 0.598 # at 298K in kcal/mol 

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


#get the actual Iteration from logger module
#iteration =
### SETUP TEST ITERATION
iteration = [iteration.Iteration(0), iteration.Iteration(1), iteration.Iteration(2)]
iteration[0].generateBin(reference_iteration_id=0, 
                 reference_bin_id=0, reference_segment_id=0, 
                 target_number_of_segments=5)
iteration[0].bins[0].generateSegment(probability=1, parent_bin_id=0, parent_segment_id=0)
iteration[0].bins[0].segments[0].parent_iteration_id=0
iteration[1].generateBin(reference_iteration_id=0, 
                 reference_bin_id=0, reference_segment_id=0, 
                 target_number_of_segments=5)
iteration[1].bins[0].generateSegment(probability=0.1, parent_bin_id=0, parent_segment_id=0)
iteration[1].bins[0].generateSegment(probability=1, parent_bin_id=0, parent_segment_id=0)
iteration[1].bins[0].segments[0].parent_iteration_id=0
iteration[1].bins[0].segments[1].parent_iteration_id=0

iteration[2].generateBin(reference_iteration_id=1, 
                 reference_bin_id=0, reference_segment_id=0, 
                 target_number_of_segments=5)
iteration[2].bins[0].generateSegment(probability=0.1, parent_bin_id=0, parent_segment_id=0)
iteration[2].bins[0].generateSegment(probability=1, parent_bin_id=0, parent_segment_id=0)
iteration[2].bins[0].segments[0].parent_iteration_id=0
iteration[2].bins[0].segments[1].parent_iteration_id=0
### 

#Calculate the coordinate values and store them together with
#the trajectory probability into coordinates 
coordinates     = numpy.zeros([0,2])
coordinates_tmp = numpy.zeros([1,2]) 
for iteration_loop in iteration:
    for bin_loop in iteration_loop.bins:
        for segment_loop in bin_loop.segments:
            coordinates_tmp[0,0] = md_module.ana_calculatePMF_getCoordinate(segment_loop, cpptraj_lines)
            coordinates_tmp[0,1] = segment_loop.getProbability()  
            coordinates          = numpy.append(coordinates, coordinates_tmp, axis=0)
 

#Calculate the weighted histogram and PMF     
#Setup variables
hist_min =  min(coordinates[:,0])
hist_max =  max(coordinates[:,0])
dcoord   =  1.0 * (hist_max - hist_min ) / args.number_of_bins
hist     =  numpy.zeros([args.number_of_bins,3], float)
#Sort coords into histogramm
for i in range(0,len(coordinates[:,0])):
    index       = int( (coordinates[i,0] - hist_min) / dcoord )
    #maximum coord entry shall not be in an extra bin:
    if index==args.number_of_bins:
        index = index - 1
    hist[index,1] = hist[index,1] + coordinates[i,1]
#Assign the bin positions and calculate free energy:
for i in range(0,args.number_of_bins):
    hist[i,0] = hist_min + i * dcoord
    if hist[i,1]>0:
        hist[i,2]  = - kT * log(hist[i,1])
    else:
        hist[i,2]  = 'NaN'
        
#Shift minimum to zero        
pmf_min = min(hist[:,2])
for i in range(0,args.number_of_bins):
    hist[i,2] = hist[i,2] - pmf_min

#Save PMF to file
header_line = 'Coordinate Value, Probability, Free energy'
numpy.savetxt(args.output_path, hist, header = header_line)


            


