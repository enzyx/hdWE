#!/usr/bin/env python2
#
# This file is part of hdWE. 
# Copyright (C) 2016 Manuel Luitz <manuel.luitz@tum.de>
# Copyright (C) 2016 Rainer Bomblies <r.bomblies@tum.de>
# Copyright (C) 2016 Fabian Zeller
#
# hdWE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# hdWE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with hdWE. If not, see <http://www.gnu.org/licenses/>.
# 
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
from lib.bin_classifier import getCoordinateIds
from lib.functions_ana_general import weightedHistogram


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
    if (segment.getIterationId(), segment.getBinId(), segment.getId()) in parent_segment_list:
        return True
    else:
        return False
    
def calculateHistogram(segment_list, number_of_histo_bins):
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
    dcoord                 = 1.0 * (hist_max - hist_min ) / number_of_histo_bins
    histbin_coordinate     = numpy.zeros([number_of_histo_bins], float)
    free_energy            = numpy.zeros([number_of_histo_bins], float)
    probability            = numpy.zeros([number_of_histo_bins], float)
    
    #Sort coords into histogram
    for data in datapoints:
        hist_index       = int( (data.coordinate - hist_min) / dcoord )
        #maximum coord entry shall not be in an extra bin:
        if hist_index==number_of_histo_bins:
            hist_index = hist_index - 1
        probability[hist_index] += data.probability
    
    #Assign the bin positions and calculate free energy:
    for i in range(0,number_of_histo_bins):
        histbin_coordinate[i] = hist_min + dcoord/2 + i * dcoord
        if probability[i]>0:        
            free_energy[i]  = - constants.kT * log(probability[i])
        else:
            free_energy[i]  = 'Inf'
            
    #Shift minimum to zero        
#     pmf_min = min(free_energy)
#     for i in range(0,number_of_histo_bins):
#         free_energy[i] -= pmf_min
    
    histogram = (histbin_coordinate, free_energy, probability)
    return histogram
    
###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    default="log", metavar="DIR", required=True,
                    help="The log directory")
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
parser.add_argument('-o', '--output', dest="outfile", 
                    type=str, default='ana_BinPMF.pmf',
                    help="Filename of output PMF data.")                          
parser.add_argument('--bin', dest="bin_id", nargs='?',
                    type=int, default=False, const=0, metavar="INT",
                    help="Bin regarded for PMF calculation.")   
parser.add_argument('-d', '--dimension', dest="dimension", default=0,
                    type=int, metavar="INT",
                    help="Dimension to calculate PMF in..")   
parser.add_argument('-N', '--number_of_histo_bins', dest="number_of_histo_bins", nargs='?',
                    default=False, type=int, const=32, 
                    help="Number of bins used to calculate the probability histogram.")
    
# Initialize
sys.stdout.write('\033[1mCalculating Bin PMF\033[0m (Free Energy is given in kcal/mol at 298K).\n')      
args = parser.parse_args()
# parse outfile input
split_outfile = args.outfile.split('.')
if len(split_outfile) > 1:
    outfileend    = split_outfile[-1]
    outfilestring = ".".join(split_outfile[:-1])
else:
    outfileend = 'pmf'
    outfilestring = args.outfile
    
if args.first_iteration - args.last_iteration == 0:
    raise Exception('Need more than one iteration since parents are compared.\n')
    

if args.plot:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
logger = Logger(args.logdir)

#get the Iterations with coordinate files from logger module
sys.stdout.write(" Loading iterations\n")
if args.last_iteration == -1:
    args.last_iteration = logger.getLastIterationId()

#Calculate the coordinate values and store them together with
#the trajectory probability into coordinates 
coords         = {'run': [], 'initial':[]} # coordinates of segments
probs          = {'run': [], 'initial':[]} # probabilities of segments
n_segments            = 0
n_initial_segments    = 0
N_ITERATIONS = args.last_iteration - args.first_iteration

        
# read in initial_segment coordinates and probabilities
prev_iteration = logger.loadIteration(args.first_iteration)
for iteration_id in range(args.first_iteration + 1, args.last_iteration+1):
    iteration = logger.loadIteration(iteration_id)
    if len(iteration.bins) <= args.bin_id:
        continue
    
    if args.bin_id:
        bins_to_scan = [iteration.bins[args.bin_id]]
    else:
        bins_to_scan = iteration.bins
    sys.stdout.write(' Reading coordinates for segments of iteration '\
                     '{it_id:08d}/{first_it:05d}-{last_it:05d}'\
                     '\r'.format(it_id     = iteration_id,
                                 first_it  = args.first_iteration,
                                 last_it   = args.last_iteration))
    sys.stdout.flush()
    
    for this_bin in bins_to_scan:
        if len(coords['run']) < this_bin.getId() + 1:
            coords['run'].append([])
            probs['run'].append([])
        for this_segment in this_bin:
            coords['run'][this_bin.getId()].append( prev_iteration.bins[this_segment.getParentBinId()]\
                                                   .segments[this_segment.getParentSegmentId()]\
                                                   .getCoordinates()\
                                                   [args.dimension] )
            probs['run'][this_bin.getId()].append(this_segment.getProbability())
            n_segments += 1
        
        if len(coords['initial']) < this_bin.getId() + 1:
            coords['initial'].append([])
            probs['initial'].append([])    
        for this_initial_segment in this_bin.initial_segments:
            coords['initial'][this_bin.getId()].append( prev_iteration.bins[this_initial_segment.getParentBinId()]\
                                                        .segments[this_initial_segment.getParentSegmentId()]\
                                                        .getCoordinates()\
                                                        [args.dimension] )
            probs['initial'][this_bin.getId()].append(this_initial_segment.getProbability())
            n_initial_segments += 1
                
    prev_iteration = iteration
    
N_BINS = len(coords['run'])
sys.stdout.write("\n")
sys.stdout.write(' found {} initial and {} run segments\n'.format(n_initial_segments, n_segments))

sys.stdout.write("\n")

# normalize probabilities:
for bin_id in range(N_BINS):
    probs['run'][bin_id]     = numpy.asarray(probs['run'][bin_id]) / N_ITERATIONS
    probs['initial'][bin_id] = numpy.asarray(probs['initial'][bin_id]) / N_ITERATIONS
    
       
#Calculate the weighted histogram and PMF
sys.stdout.write(' Creating histograms\n') 
# get number of histogram bins if not entered:
if not args.number_of_histo_bins:
    args.number_of_histo_bins = iteration.bins[0].getNumberOfSegments()
if args.plot:
    f,ax = plt.subplots(1,1)
    colors = ['Blue', 'Red']



for set_index, this_set in enumerate(['run', 'initial']):
    set_title = "{} segments".format(this_set)
    # open corresponding file
    outfilename = '{fn}.{set}.{fnend}'.format(fn = outfilestring, 
                                       set = this_set,
                                       fnend = outfileend)
    pmffile = open(outfilename, "w")
    header_line = '#Coordinate Value, Free Energy, Probability\n'
    pmffile.write(header_line)
                
    for bin_id in range(N_BINS):                 
        # create histogram
        #histogram =  weightedHistogram(zip(coords[this_set][bin_id], probs[this_set][bin_id]), args.number_of_histo_bins)
        prob_histo = numpy.histogram(coords[this_set][bin_id], 
                                    weights = probs[this_set][bin_id],
                                    bins = args.number_of_histo_bins)
        """ 
            prob_histo has as first array values of histogram
            and as second the bin edges, i.e. one more element than
            first array!
        """ 
        # get centers of bins
        bin_centers = []
        old_border = prob_histo[1][0]
        for border in prob_histo[1][1:]:
            bin_centers.append((old_border + border) / 2.0)
            
        # calculate free energies
        #  numpy matrix histogram with coordinate, free energy and probability:
        histogram = numpy.zeros([len(prob_histo[0]),3], dtype=float) 
        histogram[:,0] = bin_centers[:]
        # add (normalized) probability
        histogram[:,2] = numpy.array(prob_histo[0][:])
        for index,prob in enumerate(histogram[:,2]):
            if prob > 0:
                histogram[index,1] = -0.593 * log(prob)
            else:
                histogram[index,1] = 'Inf'
        data_to_save = histogram
        numpy.savetxt(pmffile, data_to_save)
        pmffile.write("\n")
        
        # plot directly
        if args.plot:
            #ax.scatter(histogram[0], histogram[1], color=color)
            ax.plot(histogram[:,0],
                    histogram[:,1],
                    color=colors[set_index],
                    label='{setname}'.format(setname = set_title))
    
        
    pmffile.close()
    
sys.stdout.write('\n PMF data written to: {}\n'.format(outfilename)) 

if args.plot:
    ax.legend()
    plt.show()




