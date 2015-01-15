import argparse
import numpy
import iteration
from amber_module import MD_module

parser = argparse.ArgumentParser(description=
    'Calculates the PMF along an arbitrary coordinate from the data of a hdWE run.')
parser.add_argument('-d', '--dir', type=str, 
                    dest="work_dir", required=True, metavar="DIR",
                    help="The working direcory")
parser.add_argument('-f', '--first_it', 
                    dest="first_iteration", required=False, 
                    type=int, default=0, help="First iteration to load.")                    
parser.add_argument('-l', '--last_it', 
                    dest="last_iteration", required=False, 
                    type=int, default=1, help="Last iteration to load.")  
parser.add_argument('-o', '--output',  
                    dest="output_path", required=False, 
                    type=str, default='ana_calculatePMF.output', help="output filename")  
parser.add_argument('-b', '--number_of_bins', 
                    dest="number_of_bins", required=False, 
                    type=int, default=100, help="Number of bins used to calculate the probability histogram.")  
parser.add_argument('-c', '--conf', type=str, dest="input_md_conf", 
                    required=True, metavar="FILE",
                    help="MD Software configuration file")
parser.add_argument('-i', '--cpptraj_lines_file', type=str, dest="cpptraj_lines_file_path", 
                    required=True, metavar="FILE",
                    help="cpptraj syntax line that defines the reaction coordinate.")
               
# Initialize
args = parser.parse_args()
md_module = MD_module(args.work_dir, args.input_md_conf, debug=True)
# Load cpptraj input line
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
iteration = [iteration.Iteration(0), iteration.Iteration(1)]
iteration[0].generateBin(reference_iteration_id=0, 
                 reference_bin_id=0, reference_segment_id=0, 
                 target_number_of_segments=5)
iteration[0].bins[0].generateSegment(probability=1, parent_bin_id=0, parent_segment_id=0)
iteration[0].bins[0].segments[0].parent_iteration_id=0
iteration[1].generateBin(reference_iteration_id=0, 
                 reference_bin_id=0, reference_segment_id=0, 
                 target_number_of_segments=5)
iteration[1].bins[0].generateSegment(probability=1, parent_bin_id=0, parent_segment_id=0)
iteration[1].bins[0].generateSegment(probability=1, parent_bin_id=0, parent_segment_id=0)
iteration[1].bins[0].segments[0].parent_iteration_id=0
iteration[1].bins[0].segments[1].parent_iteration_id=0
### 
 
coordinates = numpy.array([])

for iteration_loop in iteration:
    for bin_loop in iteration_loop.bins:
        for segment_loop in bin_loop.segments:
            probability      = segment_loop.getProbability()
            coordinate_value = md_module.ana_calculatePMF_getCoordinate(segment_loop, cpptraj_lines)  
            coordinates=numpy.append(coordinates, [coordinate_value, probability] )
 
print coordinates
           
#pmf = WeightedHistogram(coordinates, number_of_bins)
#pmf = -kT log (pmf)   
            


