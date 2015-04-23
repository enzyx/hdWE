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
    def __init__(self, iteration_index, bin_index, segment_index, coordinates, probability,
                 parent_iteration_id, parent_bin_id, parent_segment_id):
        self.iteration_index     = int(iteration_index)
        self.bin_index           = int(bin_index)
        self.segment_index       = int(segment_index)
        self.coordinates         = coordinates
        self.probability         = probability
        self.parent_iteration_id = parent_iteration_id
        self.parent_bin_id       = parent_bin_id
        self.parent_segment_id   = parent_segment_id 
        
    def getIterationId(self):
        return self.iteration_index
              
    def getBinId(self):
        return self.bin_index

    def getCoordinates(self):
        return self.coordinates
           
    def getNameString(self):
        """
        @return the indices in a string following the scheme iteration_bin_segment
        """
        return "{0:05d}_{1:05d}_{2:05d}".format(self.iteration_index,
                                                self.bin_index,
                                                self.segment_index)
        
    def getParentIterationId(self):
        return self.parent_iteration_id
            
    def getParentBinId(self):
        return self.parent_bin_id
            
    def getParentSegmentId(self):
        return self.parent_segment_id
    
class BinData(object):
    
    def __init__(self, bin_id):
        self.id = bin_id
        self.segments = []
     
    def getId(self):
        return self.id 
    
class IterationData(object):
    def __init__(self, iteration_id):
        self.id = iteration_id
        self.bins = []
     
    def getId(self):
        return self.id   
    
    
### methods ### 

def isParent(segment, parent_segment_list):
    if (segment_loop.getIterationId(), segment_loop.getBinId(), segment_loop.getId()) in parent_segments:
        return True
    else:
        return False

def getSegment(iterations, iteration_id, bin_id, segment_id):
    """
    goes through the iterations array and searches for the right iteration, bin and segment id
    needs - iteration with .bins list
          - bins with .segments list
          - segment with .getNameString()
    @return segment element of bins element of iterations
    """
    found_segment = False
    for iteration in iterations:
        if iteration.getId() == iteration_id:
            if len(iteration.bins) > bin_id:
                if len(iteration.bins[bin_id].segments) > segment_id:
                    parent_segment = iteration.bins[bin_id].segments[segment_id]
                    found_segment = True
        if iteration.getId() > iteration_id:
            break
    if not found_segment:
        raise Exception("failed to find segment {sn}".format(sn=segment.getNameString()))
    return parent_segment
    
def getParentSegment(segment, iterations):
    """
    goes through the iterations array and searches for the right iteration, bin and segment
    needs - iteration with .bins list
          - bins with .segments list
          - segment with .getParentIterationId()
                         .getParentBinId()
                         .getParentSegmentId()
                         .getNameString()
    @return segment element of bins element of iterations
    """
    found_parent = False
    p_iter = segment.getParentIterationId()
    p_bin  = segment.getParentBinId()
    p_seg  = segment.getParentSegmentId()
    for iteration in iterations:
        if iteration.getId() == p_iter:
            if len(iteration.bins) > p_bin:
                if len(iteration.bins[p_bin].segments) > p_seg:
                    parent_segment = iteration.bins[p_bin].segments[p_seg]
                    found_parent = True
        if iteration.getId() > p_iter:
            break
    if not found_parent:
        raise Exception("failed to find parent {i:05d}_{b:05d}_{s:05d} for segment {sn}".format(i=segment.getParentIterationId(),
                                                                                                 b=segment.getParentBinId(),
                                                                                                 s=segment.getParentSegmentId(),
                                                                                                 sn=segment.getNameString()))
    return parent_segment
    
###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    default="log", metavar="DIR", required=True,
                    help="The log directory")
parser.add_argument('-i', '--cpptraj_lines_file', dest="cpptraj_lines_file_path", 
                    type=str,   required=True,
                    help="File containig cpptraj syntax that defines the reaction coordinate.")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0,
                    help="First iteration read in. Segments of first iteration are not considered!.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Last iteration to to use for PMF calculation.")  
parser.add_argument('-o', '--output', dest="output_path", 
                    type=str, default='ana_calculatePMF.output',
                    help="Output filename")  
parser.add_argument('-p', '--plot', dest="plot", action="store_true",
                    default=False, help="plot result directly to screen.")
parser.add_argument('-c', '--conf', dest="input_md_conf", nargs='?',
                    default=False, 
                    help="Optional MD-Software configuration file")     
parser.add_argument('-t', '--threshold', dest="coordinate_threshold",
                    type=float, default=0,
                    help="Coordinate of the presumed barrier.")
parser.add_argument('-d', '--debug', dest="debug", default=False, 
                    action='store_true', help='print additional debugging info.')                             

                    
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

parent_segments = []   # a list of (iteration, bin, segment) tuples for parent segments
child_segments = []    # a list of (iteration, bin, segment) tuples for child segments
                       # where 0th parent is 0th in parent_segments list
segment_datas = []     # a list of Segment_Data structures, all segment data
references = []        # a list of Segment_Data structures, reference segments of bins
datapoints = []        # a list of Datapoint structures with coordinate and probability

first_it_id = iterations[0].getId()
last_it_id  = iterations[-1].getId()
iteration_datas = []

#if args.plot_segments:
if True:
    # find parent segments (to be saved during the coordinate calculation)
    for bin_loop in iterations[-1]:
        bin_parent_segments = []
        for segment_loop in bin_loop.initial_segments:
            parent_segments.append( (segment_loop.getParentIterationId(), segment_loop.getParentBinId(), segment_loop.getParentSegmentId()) )
            child_segments.append( (segment_loop.getIterationId(), segment_loop.getBinId(), segment_loop.getId()) )

# read in coordinates 
for iteration_loop in iterations:
    iter_data = IterationData(iteration_loop.getId())
    for bin_loop in iteration_loop:
        iter_data.bins.append(BinData(bin_loop.getId()))
        sys.stdout.write(' Calculating coordinates for iteration '\
                         '{it_id:05d}/{first_it:05d}-{last_it:05d}, '\
                         'Bin {bin_id:05d}/{bin_total:05d}\r'.format(it_id     = iteration_loop.getId(),
                                                               first_it  = first_it_id,
                                                               last_it   = last_it_id,
                                                               bin_id    = bin_loop.getId()+1,
                                                               bin_total = iteration_loop.getNumberOfBins()))
        sys.stdout.flush()
        for segment_loop in bin_loop:
            segment_coordinates = md_module.ana_calcCoordinateOfSegment(segment_loop, cpptraj_lines, use_trajectory = False)
            segment_probability = segment_loop.getProbability()
            
            iter_data.bins[-1].segments.append(SegmentData(iteration_index     = iteration_loop.getId(),
                                                           bin_index           = bin_loop.getId(),
                                                           segment_index       = segment_loop.getId(),
                                                           coordinates         = segment_coordinates[-1],
                                                           probability         = segment_probability,
                                                           parent_iteration_id = segment_loop.getParentIterationId(), 
                                                           parent_bin_id       = segment_loop.getParentBinId(), 
                                                           parent_segment_id   = segment_loop.getParentSegmentId()))
    iteration_datas.append(iter_data)
sys.stdout.write("\n")

# compare segments to parent
sys.stdout.write("iteration \t probability(events)\n")
for iteration in iterations[2:]:
    low_to_high = 0.0
    low_to_high_counter = 0
    high_to_low = 0.0
    high_to_low_counter = 0
    
    for bin_loop in iteration:
        for initial_segment in bin_loop.initial_segments:
            segment = getSegment(iteration_datas, 
                                 initial_segment.getParentIterationId(),
                                 initial_segment.getParentBinId(),
                                 initial_segment.getParentSegmentId())
            parent_segment = getParentSegment(segment, iteration_datas)
            # low to high flow
            if parent_segment.getCoordinates() < args.coordinate_threshold\
               and segment.getCoordinates() > args.coordinate_threshold:
                low_to_high += parent_segment.probability
                low_to_high_counter += 1
                if args.debug:
                    print ("{seg}({prob},{coord}) <- {pseg}({pprob},{pcoord})".format(pseg=parent_segment.getNameString(),
                                                                                 pprob=parent_segment.probability,
                                                                                 pcoord=parent_segment.getCoordinates(),
                                                                                 seg=segment.getNameString(),
                                                                                 prob=segment.probability,
                                                                                 coord=segment.getCoordinates()))
            # high to low flow
            if parent_segment.getCoordinates() > args.coordinate_threshold\
               and segment.getCoordinates() < args.coordinate_threshold:
                high_to_low += segment.probability
                high_to_low_counter += 1
                if args.debug:
                    print ("{pseg}({pprob},{pcoord}) -> {seg}({prob},{coord})".format(pseg=parent_segment.getNameString(),
                                                                                 pprob=parent_segment.probability,
                                                                                 pcoord=parent_segment.getCoordinates(),
                                                                                 seg=segment.getNameString(),
                                                                                 prob=segment.probability,
                                                                                 coord=segment.getCoordinates()))
    # iteration id is the one before because initial_segments are considered.
    # so a transition happened at the iteration before
    sys.stdout.write("{it:05d}:\t-> {plh:6f}({lhc:3d})\t|\t<- {phl:6f}({hlc:3d})\n".format(it  =iteration.getId()-1,
                                                                          plh =low_to_high,
                                                                          lhc =low_to_high_counter,
                                                                          phl =high_to_low,
                                                                          hlc =high_to_low_counter))
    #sys.stderr.write("\n")

