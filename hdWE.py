#!/usr/bin/python3

import argparse
import numpy
import initiate
from iteration import Iteration
from bin import Bin
from segment import Segment
from logger import Logger

#############################
# parsing the commandline
parser = argparse.ArgumentParser(description=
    'hdWE is a hyperdimensional weighted ensemble simulation implementation.')
parser.add_argument('-d', '--dir', type=str, 
                    dest="work_dir", required=True, metavar="DIR",
                    help="The working direcory")
#~ parser.add_argument('-t', '--top', type=str, dest="input_md_topol", 
                    #~ required=True, metavar="FILE",
                    #~ help="System topology file")
#~ parser.add_argument('-p', '--parm', type=str, dest="input_md_param", 
                    #~ required=True, metavar="FILE",
                    #~ help="MD paramter file")
parser.add_argument('-c', '--conf', type=str, dest="input_md_conf", 
                    required=True, metavar="FILE",
                    help="The starting structure file")
parser.add_argument('--segments-per-bin', type=int, dest="segments_per_bin", 
                    metavar="10", default=10, nargs='?',
                    help="Number of trajectories per bin")
parser.add_argument('--iterations', type=int, dest="max_iterations", 
                    metavar="50", default=50, nargs='?',
                    help="Maximum number of iterations")
parser.add_argument('--threshold', type=float, dest="coordinate_threshold", 
                    metavar="0.1", default=0.1, nargs='?',
                    help="Defines the minimal RMSD of a trajectory to all other bins "
                         "after which a new bin is created")
parser.add_argument('--minimal-probability', type=float, dest="input_minimal_probability", 
                    metavar="0.01", default=0.01, nargs='?',
                    help="Minimal probability a trajectory must have to"
                    " allow forking a new bin")
args = parser.parse_args()
#############################

#############################
# The global list of arrays
iterations = []

# Initialize the logger
logger = Logger("logfile.log")

initiate.prepare(args.work_dir, starting_sturcture="", override="", debug="")
iterations.append(initiate.create_initial_iteration(args.segments_per_bin))

#~ if("amber"):
    #~ if self.MD_software=='amber':
from amber_module import MD_module

md_module = MD_module(args.work_dir, args.input_md_conf, debug=False)

#~ logger.read(logfile)

# Loop
for iteration_counter in range(1, args.max_iterations):
    iteration = Iteration(iteration_counter)
    parent_iteration = iterations[iteration_counter - 1]
    # Generate all previous bins for new iteration
    for parent_bin in parent_iteration:
        iteration.generateBin(reference_iteration_id=parent_bin.getReferenceIterationId(),
                              reference_bin_id=parent_bin.getReferenceBinId(),
                              reference_segment_id=parent_bin.getReferenceSegmentId(),
                              target_number_of_segments=args.segments_per_bin)
    
    coordinates = numpy.array([]) # numpy array?
    # Sort segments in bins. Generate new bin if required
    for parent_bin in parent_iteration:
        for segment in parent_bin:
            coordinates = md_module.CalculateCoordinate(segment, iteration.bins)
            min_coordinate = numpy.min(coordinates)
            # Sort Segment into appropriate bin
            if (min_coordinate <= args.coordinate_threshold):
                bin_id = numpy.argmin(coordinates)
                iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                  parent_bin_id=segment.getBinId(),
                                                  parent_segment_id=segment.getId())
            # If necessary create new bin
            else:
                bin_id = iteration.generateBin(reference_iteration_id=segment.getIterationId(),
                                      reference_bin_id=segment.getBinId(),
                                      reference_segment_id=segment.getId(),
                                      target_number_of_segments=args.segments_per_bin)
                iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                  parent_bin_id=segment.getBinId(),
                                                  parent_segment_id=segment.getId())
                                                  
    # Split and merge (Manu)
    for this_bin in iteration:
        this_bin.resampleSegments()
                
    # log iteration (Rainer)
    logger.log_iterations(iterations[-1:])
    
    
    
    # Run MD
    md_module.RunMDs(iteration)
    

                                    
    iterations.append(iteration)  
    
logger.close()
