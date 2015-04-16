"""
Calculates the PMF along an arbitrary coordinate 
from the data of a hdWE run.
"""
from __future__ import print_function
import numpy
import matplotlib.pyplot as plt  
import constants
import sys
import segment
from logger import Logger
from math import log
from amber_module import MD_module
import argparse  

### classes ### 
class SegmentData(object):
    """
    a container for resulting data of a segment
    """
    def __init__(self, iteration_index, bin_index, coordinate):
        self.iteration_index = iteration_index
        self.bin_index = bin_index
        self.coordinate = coordinate
        
    def getIterationId(self):
        return self.iteration_index
              
    def getBinId(self):
        return self.bin_index
              
    def getCoordinate(self):
        return self.coordinate

###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-c', '--conf', dest="input_md_conf", 
                    type=str, required=True,
                    help="MD-Software configuration file")
parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    default="log", metavar="DIR",
                    help="The log directory")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0,
                    help="First iteration to use for PMF calculation.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Last iteration to to use for PMF calculation.")  
parser.add_argument('-o', '--output', dest="output_path", 
                    type=str, default='ana_calculatePMF.output',
                    help="Output filename")  
parser.add_argument('--op', '--plotfile', dest="output_plot", 
                    type=str, default='ana_calculatePMF.pdf',
                    help="Filename for output plot.")  
parser.add_argument('-N', '--number_of_bins', dest="number_of_bins",
                    type=int, default=100, 
                    help="Number of bins used to calculate the probability histogram.")  
parser.add_argument('-i', '--cpptraj_lines_file', dest="cpptraj_lines_file_path", 
                    type=str,   required=True,
                    help="File containig cpptraj syntax that defines the reaction coordinate.")
parser.add_argument('-p', '--plot', dest="plot", action="store_true",
                    default=False, help="plot result directly to screen.")    
parser.add_argument('-r', '--reference', nargs='?', default=False, 
                    dest="reference", metavar="FILE", const="reference.dat",
                    help="The (optional) reference file for comparison with cMD.")  
parser.add_argument('--rc', '--reference_column', nargs='?', default=False, 
                    dest="ref_column", metavar="INT", const=2, type=int,
                    help="The (optional) reference file for comparison with cMD.")                      
parser.add_argument('-y', dest="use_trajectory", action="store_true",
                    default=False, help="use conformations from the trajectories,"\
                                        "not only the end conformation")                                                   
parser.add_argument('--segments', dest="plot_segments", action="store_true",
                    default=False, help="plot final segment positions of last iteration.")
parser.add_argument('--bin-references', dest="plot_bin_references", action="store_true",
                    default=False, help="plot bin reference positions")                                        
                    
# Initialize
print('\033[1mCalculating PMF\033[0m (Free Energy is given in kcal/mol at 298K).')      
args = parser.parse_args()
md_module = MD_module(args.input_md_conf, debug=False)

#get the actual Iteration from logger module
logger = Logger(args.logdir)
iterations = logger.loadIterations(begin=args.first_iteration, end=args.last_iteration)

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

segments = []
references = []
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
            coordinate = md_module.ana_calcCoordinateOfSegment(segment_loop, cpptraj_lines, use_trajectory = args.use_trajectory)
            if iteration_loop.getId() == iterations[-1].getId():
                segments.append(SegmentData(iteration_loop.getId(), bin_loop.getId(), coordinate))
            coordinates_tmp_tmp = coordinate
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
    coordinate = md_module.ana_calcCoordinateOfSegment(segment_tmp, cpptraj_lines, False)
    references.append(SegmentData(iteration_loop.getId(), bin_loop.getId(), coordinate))
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
print('\n PMF data written to: ' + args.output_path) 

# Plotting
segment_colors = ["Blue", "Red", "Green", "Orange", "c", "m", "y"]
seg_step = 0.3
f, ax = plt.subplots(1,1)
ax.grid()

# load reference (cMD) data
if args.reference:
    ref_data = numpy.transpose(numpy.loadtxt(args.reference))
    ax.plot(ref_data[0], ref_data[args.ref_column-1], label=args.reference)

ax.plot(hist[:,0], hist[:,1], label=args.logdir)

if args.plot_segments:
    for bin_id in range(iterations[-1].getNumberOfBins()):
        seg_x = []
        seg_y = []
        for segment_loop in segments:
            if segment_loop.getBinId() == bin_id:
                seg_x.append(segment_loop.getCoordinate())
                seg_y.append(-seg_step - seg_step*segment_loop.getBinId())
                
            if segment.getBinId() > bin_id:
                break
        ax.scatter(seg_x, 
                seg_y,
                marker="s",
                color=segment_colors[bin_id%len(segment_colors)])

if args.plot_bin_references:
    bin_ref_x = []
    bin_ref_y = []    
    for index, bin_loop in enumerate(iterations[-1]):
        bin_ref_x.append(references[index].getCoordinate())
        bin_ref_y.append(-1)
        bin_color = 1.0*bin_loop.getId()/iterations[-1].getNumberOfBins()
    cmap = plt.get_cmap("coolwarm")
    ax.scatter(bin_ref_x, 
               bin_ref_y,
               marker="s")               

ax.legend()
if args.plot:
    plt.show() 
else:
    plt.savefig(args.output_plot, bbox_inches='tight', transparent=True)
    print('\n Output written to: ' + args.output_plot)

   

