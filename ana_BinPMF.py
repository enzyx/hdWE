#!/usr/bin/env python2
"""
Calculates the PMF along an arbitrary coordinate 
for initial segments and actually run segments. 
"""
from __future__ import print_function
import sys
import numpy
from math import log
import lib.constants as constants
import lib.segment
from lib.logger import Logger
from lib.amber_module import MD_module
import argparse  
import ConfigParser


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
    def __init__(self, iteration_index, bin_index, coordinates, probability):
        self.iteration_index = iteration_index
        self.bin_index = bin_index
        self.coordinates = coordinates
        self.probability = probability
        
    def getIterationId(self):
        return self.iteration_index
              
    def getBinId(self):
        return self.bin_index

    def getCoordinates(self):
        return self.coordinates
    
    def getProbability(self):
        return self.probability
    
### methods ### 

def isParent(segment, parent_segment_list):
    if (segment_loop.getIterationId(), segment_loop.getBinId(), segment_loop.getId()) in parent_segments:
        return True
    else:
        return False
    
def calculateHistogram(segment_list, number_of_bins):
    """
    Calculates a histogram of the SegmentData objects with a probability and a coordinates list
    """
    datapoints = []
    for segment_data in segment_list:
        for coordinate in segment_data.getCoordinates():
            datapoints.append(Datapoint(coordinate, 
                                        segment_data.getProbability()/len(segment_data.getCoordinates())))
    datapoints             = sorted(datapoints, key=lambda datapoint: datapoint.coordinate)
    hist_min               = datapoints[0].coordinate
    hist_max               = datapoints[-1].coordinate
    dcoord                 = 1.0 * (hist_max - hist_min ) / number_of_bins
    histbin_coordinate     = numpy.zeros([number_of_bins], float)
    free_energy            = numpy.zeros([number_of_bins], float)
    probability            = numpy.zeros([number_of_bins], float)
    
    #Sort coords into histogram
    for data in datapoints:
        hist_index       = int( (data.coordinate - hist_min) / dcoord )
        #maximum coord entry shall not be in an extra bin:
        if hist_index==number_of_bins:
            hist_index = hist_index - 1
        probability[hist_index] += data.probability
    
    #Assign the bin positions and calculate free energy:
    for i in range(0,number_of_bins):
        histbin_coordinate[i] = hist_min + dcoord/2 + i * dcoord
        if probability[i]>0:        
            free_energy[i]  = - constants.kT * log(probability[i])
        else:
            free_energy[i]  = 'Inf'
            
    #Shift minimum to zero        
#     pmf_min = min(free_energy)
#     for i in range(0,number_of_bins):
#         free_energy[i] -= pmf_min
    
    histogram = (histbin_coordinate, free_energy, probability)
    return histogram
    
###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    default="log", metavar="DIR", required=True,
                    help="The log directory")
parser.add_argument('-i', '--cpptraj_lines_file', dest="cpptraj_lines_file_path", 
                    type=str,   required=True,
                    help="File containig cpptraj syntax that defines the reaction coordinate.")
parser.add_argument('-c', '--conf', dest="configfile", nargs='?',
                    default=False, 
                    help="Optional MD-Software configuration file")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0,
                    help="First iteration to use for PMF calculation.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Last iteration to use for PMF calculation.") 
parser.add_argument('-p', '--plot', dest="plot", action="store_true",
                    default=False, help="plot result directly to screen.")
parser.add_argument('-o', '--output', dest="output_path", 
                    type=str, default='ana_BinPMF.pmf',
                    help="Filename of output PMF data.")                          
parser.add_argument('--bin', dest="bin_id", 
                    type=int, default=0, metavar="INT",
                    help="Bin regarded for PMF calculation.")   
parser.add_argument('-y', dest="use_trajectory", action="store_true",
                    default=False, help="use conformations from the trajectories,"\
                                        "not only the end conformation")
parser.add_argument('-N', '--number_of_bins', dest="number_of_bins", nargs='?',
                    default=False, type=int, const=32, 
                    help="Number of bins used to calculate the probability histogram.")
parser.add_argument('-a', '--average', dest="averaged_pmfs", nargs='?',
                    default=False, type=int, const=4, 
                    help="Calculate averages and produce given number of PMFs.")
    
# Initialize
sys.stdout.write('\033[1mCalculating Bin PMF\033[0m (Free Energy is given in kcal/mol at 298K).\n')      
args = parser.parse_args()
if args.plot:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
logger = Logger(args.logdir)
keep_coords_frequency = int(logger.loadConfigParameter('keep-coords-frequency', 
                                                       iteration_id = args.first_iteration))

#get the Iterations with coordinate files from logger module
if args.last_iteration == -1:
    args.last_iteration = logger.getLastIterationId()
iterations_with_coords = [i for i in range(args.first_iteration, args.last_iteration+1) \
                                if i % keep_coords_frequency == 1]
iterations = logger.loadIterationList(iterations_with_coords)

# load md module
if not args.configfile:
    args.configfile = logger.loadConfigFile(iterations[0].getId())
md_module = MD_module(args.configfile, debug=False)

# Load cpptraj input file as one string with linebreaks and delete the last line break
try:
    cpptraj_lines_file=open(args.cpptraj_lines_file_path, 'r')
except:
    sys.stderr.write('Error: could not open ' + args.cpptraj_lines_file_path)
cpptraj_lines=''
for line in cpptraj_lines_file:
    cpptraj_lines = cpptraj_lines + line
cpptraj_lines = cpptraj_lines[0:-1]
cpptraj_lines_file.close()

#Calculate the coordinate values and store them together with
#the trajectory probability into coordinates 
segment_datas         = [] # a list of lists of SegmentData structures of length len(iterations)
initial_segment_datas = [] # the same for initial segments
n_segments            = 0
n_initial_segments    = 0

first_it_id = iterations[0].getId()
last_it_id  = iterations[-1].getId()
        
# read in initial_segment coordinates and probabilities
for iteration in iterations:
    if len(iteration.bins) <= args.bin_id:
        continue
    if iteration.getId() % keep_coords_frequency != 1:
        continue
    
    initial_segment_datas.append([]) 
    this_bin = iteration.bins[args.bin_id]
    sys.stdout.write(' Calculating coordinates for initial segments of iteration '\
                     '{it_id:05d}/{first_it:05d}-{last_it:05d}'\
                     '\r'.format(it_id     = iteration.getId(),
                                 first_it  = first_it_id,
                                 last_it   = last_it_id))
    sys.stdout.flush()
    
    for this_initial_segment in this_bin.initial_segments:
        initial_segment_coordinates = md_module.ana_calcCoordinateOfSegment(this_initial_segment.getParentNameString(), 
                                                                            cpptraj_lines,
                                                                            use_trajectory = args.use_trajectory)
        initial_segment_probability = this_initial_segment.getProbability()
        initial_segment_datas[-1].append(SegmentData(iteration.getId(),
                                                     args.bin_id,
                                                     initial_segment_coordinates,
                                                     initial_segment_probability))
        n_initial_segments += 1
sys.stdout.write("\n")
sys.stdout.write(' found {} initial segments\n'.format(n_initial_segments))


# read in coordinates and probability 
for iteration in iterations:
    if len(iteration.bins) <= args.bin_id:
        continue
    # skip iterations where the parent iteration's coordinates were not kept
    if iteration.getId() % keep_coords_frequency != 1:
        continue

    segment_datas.append([])

    this_bin = iteration.bins[args.bin_id]
    sys.stdout.write(' Calculating coordinates for run segments of iteration '\
                     '{it_id:05d}/{first_it:05d}-{last_it:05d}'\
                     '\r'.format(it_id     = iteration.getId() - 1,
                                 first_it  = first_it_id,
                                 last_it   = last_it_id))
    sys.stdout.flush()

    for this_segment in this_bin:
        segment_coordinates = md_module.ana_calcCoordinateOfSegment(this_segment.getParentNameString(), 
                                                                    cpptraj_lines,
                                                                    use_trajectory = args.use_trajectory)
        segment_probability = this_segment.getProbability()
        segment_datas[-1].append(SegmentData(iteration.getId(),
                                         args.bin_id,
                                         segment_coordinates,
                                         segment_probability))
        n_segments += 1
sys.stdout.write("\n")
sys.stdout.write(' found {} run segments\n'.format(n_segments))


#segment_datas = initial_segment_datas
sys.stdout.write("\n")
       
#Calculate the weighted histogram and PMF
sys.stdout.write(' Creating histograms\n') 
# get number of bins if not entered:
if not args.number_of_bins:
    args.number_of_bins = iterations[-1].bins[0].getNumberOfSegments() 

if args.plot:
    f,ax = plt.subplots(1,1)
    cm = plt.get_cmap('gist_rainbow')
    cNorm  = mpl.colors.Normalize(vmin=0, vmax=len(segment_datas))
    scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=cm)

if not args.averaged_pmfs:
    for set_index, data_sets in enumerate([segment_datas, initial_segment_datas]):
        if set_index == 0:
            set_title = 'run segments'
        if set_index == 1: 
            set_title = 'initial segments'
        
        # open corresponding file
        pmffile = open('{fn}.{set}'.format(fn = args.output_path,
                                           set = set_title.split()[0]), "w")
        header_line = '#Coordinate Value, Free Energy, Probability\n'
        pmffile.write(header_line)
        
        # calculate histogram
        for i,segment_list in enumerate(data_sets):
            if len(segment_list) > 0:
                iteration_id = segment_list[0].getIterationId()
                histogram =  calculateHistogram(segment_list, args.number_of_bins)
                data_to_save = numpy.transpose([histogram[0],
                                                histogram[1],
                                                histogram[2]])
                numpy.savetxt(pmffile, data_to_save)
                pmffile.write("\n")
            
            # plot directly
            if args.plot:
                color = cm(1.0*i/len(data_sets))
                color = scalarMap.to_rgba(i)
                #ax.scatter(histogram[0], histogram[1], color=color)
                ax.plot(histogram[0], 
                        histogram[1], 
                        color=color, 
                        label='{} - iteration {}'.format(set_title, iteration_id))
                
        pmffile.close()
    
if args.averaged_pmfs:
    iterations_per_pmf = int(len(iterations)/args.averaged_pmfs)
    if iterations_per_pmf == 0:
        sys.stderr.write("WARNING: Too few iterations for desired number of average pmfs.\n"\
                         "         Averaging over all iterations.\n")
        iterations_per_pmf = len(iterations)
    
    for set_index, data_sets in enumerate([segment_datas, initial_segment_datas]):
        if set_index == 0:
            set_title = 'run segments'
        if set_index == 1: 
            set_title = 'initial segments'

        # open corresponding file
        pmffile = open('{fn}.{set}'.format(fn = args.output_path, 
                                           set = set_title.split()[0]), "w")
        header_line = '#Coordinate Value, Free Energy, Probability\n'
        pmffile.write(header_line)            
            
        segments_for_average = []
        for i, segment_list in enumerate(data_sets):
            if i % iterations_per_pmf == 0 and len(iterations)-i >= iterations_per_pmf:
                segments_for_average.append([])
            segments_for_average[-1] += segment_list
                
        for i,segment_list in enumerate(segments_for_average):
            # create and save histogram
            histogram =  calculateHistogram(segment_list, args.number_of_bins)
            data_to_save = numpy.transpose([histogram[0],
                                        histogram[1],
                                        histogram[2]])
            numpy.savetxt(pmffile, data_to_save)
            pmffile.write("\n")
            
            # plot directly
            if args.plot:
                cNorm  = mpl.colors.Normalize(vmin=0, vmax=len(segments_for_average))
                scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=cm)
                color = cm(1.0*i/len(segments_for_average))
                color = scalarMap.to_rgba(i)
                #ax.scatter(histogram[0], histogram[1], color=color)
                ax.plot(histogram[0],
                        histogram[1],
                        color=color,
                        label='{setname} - <it{first:05d}-it{last:05d}>'.format(setname = set_title,
                                                                                first = segments_for_average[i][0].getIterationId(),
                                                                                last  = segments_for_average[i][-1].getIterationId()))
        
                
    pmffile.close()
    
sys.stdout.write('\n PMF data written to: {}*\n'.format(args.output_path)) 

if args.plot:
    ax.legend()
    plt.show()




