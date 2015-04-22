#!/usr/bin/env python2
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

### classes ### 
class Datapoint(object):
    """
    a container for PMF datapoints with probability and coordinate
    """
    def __init__(self, coordinate, probability):
        self.coordinate = coordinate
        self.probability = probability

class SegmentData(object):
    """
    a container for resulting data of a segment
    """
    def __init__(self, iteration_index, bin_index, coordinates):
        self.iteration_index = iteration_index
        self.bin_index = bin_index
        self.coordinates = coordinates
        
    def getIterationId(self):
        return self.iteration_index
              
    def getBinId(self):
        return self.bin_index

    def getCoordinates(self):
        return self.coordinates
    
###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-c', '--conf', dest="input_md_conf", nargs='?',
                    default=False, 
                    help="Optional MD-Software configuration file")
parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    default="log", metavar="DIR", required=True,
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
if args.plot:
    import matplotlib.pyplot as plt  

# failproofing for segment plotting
if args.first_iteration == args.last_iteration and \
   args.plot_segments:
    raise Exception("Need more than 1 iteration for --segments\n")

#get the actual Iteration from logger module
logger = Logger(args.logdir)
iterations = logger.loadIterations(begin=args.first_iteration, end=args.last_iteration)

# load md module
if not args.input_md_conf:
    args.input_md_conf = logger.loadConfigFile(iterations[0].getId())
md_module = MD_module(args.input_md_conf, debug=False)

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

parent_segments = []   # a list of (bin, segment) tuples for parent segments
segment_datas = []     # a list of Segment_Data structures, all segment data
references = []        # a list of Segment_Data structures, reference segments of bins
datapoints = []        # a list of Datapoint structures with coordinate and probability

first_it_id = iterations[0].getId()
last_it_id  = iterations[-1].getId()
data_per_segment = len(md_module.ana_calcCoordinateOfSegment(iterations[0].bins[0].segments[0],
                                                             cpptraj_lines,
                                                             use_trajectory = args.use_trajectory))

if args.plot_segments:
    # find parent segments (to be saved during the coordinate calculation)
    for bin_loop in iterations[-1]:
        for segment_loop in bin_loop:
            parent_segments.append( (segment_loop.getParentBinId(), segment_loop.getParentSegmentId()) )
    parent_segments = set(parent_segments)       

# read in coordinates and probability 
for iteration_loop in iterations:
    for bin_loop in iteration_loop:
        sys.stdout.write(' Calculating coordinates for iteration '\
                         '{it_id:05d}/{first_it:05d}-{last_it:05d}, '\
                         'Bin {bin_id:05d}/{bin_total:05d}\r'.format(it_id     = iteration_loop.getId(),
                                                               first_it  = first_it_id,
                                                               last_it   = last_it_id,
                                                               bin_id    = bin_loop.getId()+1,
                                                               bin_total = iteration_loop.getNumberOfBins()))
        sys.stdout.flush()
        for segment_loop in bin_loop:
            segment_coordinates = md_module.ana_calcCoordinateOfSegment(segment_loop, cpptraj_lines, use_trajectory = args.use_trajectory)
            segment_probability = segment_loop.getProbability()
            for i,coord in enumerate(segment_coordinates):
                datapoints.append(Datapoint(coord, segment_probability))
                
            # store parent segment data
            if args.plot_segments and iteration_loop.getId() == last_it_id - 1:
                segment_datas.append(SegmentData(segment_loop.getIterationId(), bin_loop.getId(), segment_coordinates))
sys.stdout.write("\n")

# Calculate the coordinate values of the bin reference structures
sys.stdout.write(' Calculating coordinates of bin references\n')
for bin_loop in iterations[-1]:
    #create temporary segment to pass to md_module, because the bin reference segment
    #is not necessarilly within the loaded range of iterations
    segment_tmp = segment.Segment(probability = 0, parent_bin_id = 0, parent_segment_id = 0,
                  iteration_id = bin_loop.getReferenceIterationId(),
                  bin_id       = bin_loop.getReferenceBinId(),
                  segment_id   = bin_loop.getReferenceSegmentId() )
    coordinate = md_module.ana_calcCoordinateOfSegment(segment_tmp, 
                                                       cpptraj_lines, 
                                                       use_trajectory=False)
    references.append(SegmentData(iteration_loop.getId(), bin_loop.getId(), coordinate))       
            
#Calculate the weighted histogram and PMF 
sys.stdout.write(' Creating weighted histogram and PMF\n')
#Setup variables
datapoints                  = sorted(datapoints, key=lambda datapoint: datapoint.coordinate)
hist_min                    = datapoints[0].coordinate
hist_max                    = datapoints[-1].coordinate
dcoord                      = 1.0 * (hist_max - hist_min ) / args.number_of_bins
hist_histbin_coordinate     = numpy.zeros([args.number_of_bins], float)
hist_free_energy            = numpy.zeros([args.number_of_bins], float)
hist_probability            = numpy.zeros([args.number_of_bins], float)
hist_datapoints_per_histbin = numpy.zeros([args.number_of_bins], int)
hist_references_per_histbin = numpy.zeros([args.number_of_bins], int)

#Sort coords into histogram
for data in datapoints:
    hist_index       = int( (data.coordinate - hist_min) / dcoord )
    #maximum coord entry shall not be in an extra bin:
    if hist_index==args.number_of_bins:
        hist_index = hist_index - 1
    hist_probability[hist_index] += data.probability/data_per_segment
    hist_datapoints_per_histbin[hist_index] += 1

#Sort bin coords into histogram
for reference in references:
    hist_index = int( (reference.getCoordinates()[-1] - hist_min) / dcoord )
    #maximum bin coord entry shall not be in an extra bin:
    if hist_index>=args.number_of_bins:
        hist_index = args.number_of_bins - 1
    hist_references_per_histbin[hist_index] += 1
#Assign the bin positions and calculate free energy:
for i in range(0,args.number_of_bins):
    hist_histbin_coordinate[i] = hist_min + dcoord/2 + i * dcoord
    if hist_probability[i]>0:
        hist_free_energy[i]  = - constants.kT * log(hist_probability[i])
    else:
        hist_free_energy[i]  = 'Inf'
        
#Shift minimum to zero        
pmf_min = min(hist_free_energy)
for i in range(0,args.number_of_bins):
    hist_free_energy[i] -= pmf_min

#Save PMF to file
header_line = 'Coordinate Value, Free Energy, Probability, Values per hist-Bin, hdWE Bins per hist-Bin'
data_to_save = numpy.transpose([hist_histbin_coordinate,
                                hist_free_energy,
                                hist_probability,
                                hist_datapoints_per_histbin,
                                hist_references_per_histbin])
numpy.savetxt(args.output_path, data_to_save, header = header_line)
print('\n PMF data written to: ' + args.output_path) 


# Plotting
if args.plot:
    sys.stdout.write(' Preparing plot\n')
    try:
        segment_colors = ["Blue", "Red", "Green", "Orange", "c", "m", "y"]
        seg_step = 0.3
        f, ax = plt.subplots(1,1)
        ax.grid()
        
        # plot reference (cMD) data
        if args.reference:
            ref_data = numpy.transpose(numpy.loadtxt(args.reference))
            ax.plot(ref_data[0], ref_data[args.ref_column-1], label=args.reference)
        
        # plot calculated PMF    
        ax.plot(hist_histbin_coordinate, hist_free_energy, label=args.logdir)
        
        if args.plot_segments:
            for bin_id in range(iterations[-1].getNumberOfBins()):
                seg_x = []
                seg_y = []
                for segment_data in segment_datas:
                    if segment_data.getBinId() == bin_id:
                        # use last coordinate since that should be the sorted structure
                        seg_x.append(segment_data.getCoordinates()[-1])
                        seg_y.append(-seg_step - seg_step*segment_data.getBinId())
                      
                    if segment_data.getBinId() > bin_id:
                        break
                ax.scatter(seg_x, 
                           seg_y,
                           marker="s",
                           color=segment_colors[bin_id%len(segment_colors)])
        
        if args.plot_bin_references:
            bin_ref_x = []
            bin_ref_y = []    
            for hist_index, bin_loop in enumerate(iterations[-1]):
                bin_ref_x.append(references[hist_index].getCoordinates()[-1])
                bin_ref_y.append(-1)
                bin_color = 1.0*bin_loop.getId()/iterations[-1].getNumberOfBins()
            cmap = plt.get_cmap("coolwarm")
            ax.scatter(bin_ref_x, 
                    bin_ref_y,
                    marker="s")               
        
        ax.legend()
        
        plt.savefig(args.output_plot)
        print(' Plot written to to:  ' + args.output_plot)
        plt.show() 
    except:
        sys.stderr.write(" Plotting with mathplotlib failed.")

   

