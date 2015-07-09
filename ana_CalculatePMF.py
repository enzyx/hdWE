#!/usr/bin/env python2
"""
Calculates the PMF along an arbitrary coordinate 
from the data of a hdWE run.
"""
from __future__ import print_function
import numpy
import lib.constants as constants
import sys
from lib.logger import Logger
from math import log
from lib.amber_module import MD_module
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
    
### methods ### 

def isParent(segment, parent_segment_list):
    if (segment_loop.getIterationId(), segment_loop.getBinId(), segment_loop.getId()) in parent_segments:
        return True
    else:
        return False
    
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
                    default=False, help="plot last binned segment positions \n.")
parser.add_argument('--bin-references', dest="plot_bin_references", action="store_true",
                    default=False, help="plot bin reference positions")                                        

                    
# Initialize
print('\033[1mCalculating PMF\033[0m (Free Energy is given in kcal/mol at 298K).')      
args = parser.parse_args()
first_it_id = args.first_iteration
last_it_id  = args.last_iteration
if args.plot:
    import matplotlib.pyplot as plt  

# failproofing for segment plotting
if first_it_id == last_it_id and \
   args.plot_segments:
    raise Exception("Need more than 1 iteration for --segments\n")

#get the actual Iteration from logger module
logger = Logger(args.logdir)
keep_coords_frequency = int(logger.loadConfigParameter('keep-coords-frequency', 
                                                       iteration_id = first_it_id))
first_iteration = logger.loadIteration(first_it_id)
if last_it_id == -1:
    last_it_id= logger.getLastIterationId()

# load md module
if not args.input_md_conf:
    args.input_md_conf = logger.loadConfigFile(first_iteration.getId())
md_module = MD_module(args.input_md_conf, debug=False)

#Calculate the coordinate values and store them together with
#the trajectory probability into coordinates 

parent_segments = []   # a list of (iteration, bin, segment) tuples for parent segments
child_segments = []    # a list of (iteration, bin, segment) tuples for child segments
                       # where 0th parent is 0th in parent_segments list
segment_datas = []     # a list of Segment_Data structures, all segment data
references = []        # a list of Segment_Data structures, reference segments of bins
datapoints = []        # a list of Datapoint structures with coordinate and probability

if args.plot_segments:
    # find parent segments (to be saved during the coordinate calculation)
    # of last iteration where the parent coordinates are stored
    segments_iteration_id = (last_it_id-1) - ((last_it_id-1) % keep_coords_frequency)
    segments_iteration = logger.loadIteration(segments_iteration_id + 1)
    for bin_loop in segments_iteration:
        bin_parent_segments = []
        for segment_loop in bin_loop.initial_segments:
            parent_segments.append( (segment_loop.getParentIterationId(), segment_loop.getParentBinId(), segment_loop.getParentSegmentId()) )
            child_segments.append( (segment_loop.getIterationId(), segment_loop.getBinId(), segment_loop.getId()) )

# read in coordinates and probability
for i in range(first_it_id, last_it_id + 1):
    # skip iterations without saved data
    if i%keep_coords_frequency != 0:
        continue
    
    sys.stdout.write(' Reading coordinates for iteration '\
                     '{it_id:05d}/{first_it:05d}-{last_it:05d}\r'.format(it_id     = i,
                                                           first_it  = first_it_id,
                                                           last_it   = last_it_id))
    sys.stdout.flush()
    current_iteration = logger.loadIteration(i)
        
    for bin_loop in current_iteration:
        for segment_loop in bin_loop:
            segment_coordinates = segment_loop.getCoordinates()
            segment_probability = segment_loop.getProbability()
            for i,coord in enumerate(segment_coordinates):
                datapoints.append(Datapoint(coord, segment_probability))
                
            # store parent segment data
            if args.plot_segments:
                if isParent(segment_loop, parent_segments):
                    index = parent_segments.index((segment_loop.getIterationId(), segment_loop.getBinId(), segment_loop.getId()) )
                    segment_datas.append(SegmentData(child_segments[index][0], child_segments[index][1], segment_coordinates))
sys.stdout.write("\n")       
last_iteration = current_iteration
# the number of datapoints per segment
# (if trajectory was used and more than one configuration is available)                    
data_per_segment = len(segment_coordinates)                     
                
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
    #try:
    if True:    
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
            for bin_id in range(last_iteration.getNumberOfBins()):
                seg_x = []
                seg_y = []
                for segment_data in segment_datas:
                    if segment_data.getBinId() == bin_id:
                        # use last coordinate since that should be the sorted structure
                        seg_x.append(segment_data.getCoordinates()[-1])
                        seg_y.append(-seg_step - seg_step*segment_data.getBinId())                  

                ax.scatter(seg_x, 
                           seg_y,
                           marker="s",
                           color=segment_colors[bin_id%len(segment_colors)])                 
        
        if args.plot_bin_references:
            bin_ref_x = []
            bin_ref_y = []    
            for hist_index, bin_loop in enumerate(last_iteration):
                bin_ref_x.append(references[hist_index].getCoordinates()[-1])
                bin_ref_y.append(-1)
                bin_color = 1.0*bin_loop.getId()/last_iteration.getNumberOfBins()
            cmap = plt.get_cmap("coolwarm")
            ax.scatter(bin_ref_x, 
                    bin_ref_y,
                    marker="s")               
        
        ax.legend()
        
        plt.savefig(args.output_plot)
        print(' Plot written to to:  ' + args.output_plot)
        plt.show() 
    #except:
        #sys.stderr.write(" Plotting with matplotlib failed.")

   

