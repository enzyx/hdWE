#!/usr/bin/python3

import sys
import os
import argparse
import numpy
import initiate
import analysis_operations
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
parser.add_argument('-c', '--conf', type=str, dest="input_md_conf", 
                    required=True, metavar="FILE",
                    help="The starting structure file")
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")                     
parser.add_argument('--segments-per-bin', type=int, dest="segments_per_bin", 
                    metavar="50", default=10, nargs='?',
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
parser.add_argument('--debug', dest="debug", action="store_true",
                    default=False, help="Turn debugging on")
parser_mdgroup = parser.add_mutually_exclusive_group(required=True)
parser_mdgroup.add_argument("--amber", dest="amber", action="store_true",
                    default=False)
parser_mdgroup.add_argument("--gromacs", dest="gromacs", action="store_true",
                    default=False)
args = parser.parse_args()
# guarantee a working work_dir variable
if args.work_dir[-1] != "/":
    args.work_dir +="/"
print (args.work_dir)
#############################
    
# The global list of arrays
iterations = [] 

# Initialize the logger

logger = Logger(args.work_dir + args.logfile)
# Read previous log
#~ print(os.stat(args.logfile).st_size)
#~ if os.stat(args.logfile).st_size != 31133:
    #~ iterations = logger.load_iterations()

# Setup the work_dir and initiate iterations
initiate.prepare(args.work_dir, starting_structure="", override="", debug=args.debug)
if len(iterations)==0:
    iterations.append(initiate.create_initial_iteration(args.segments_per_bin))

# Check MD suite
if(args.amber):
    from amber_module import MD_module
    md_module = MD_module(args.work_dir, args.input_md_conf, debug=args.debug)
if(args.gromacs):
    print("Sorry, support for gromacs is not implemented yet.")
    sys.exit(-1)


# Loop
for iteration_counter in range(len(iterations), args.max_iterations):
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
                
    if args.debug: 
        print("The overall probabiliy is at {0:05f}".format(iteration.getProbability()))
    
    # Run MD
    md_module.RunMDs(iteration)
    
    iterations.append(iteration)
    
    #test: print rate matrix, bin probabilities from rates and actual bin probabilities
    if iteration_counter > 15:
        mean_rate_matrix= analysis_operations.meanRateMatrix(iterations, iteration_counter -10, iteration_counter)
        print('')
        try:
            print('Bin probabilities from Rate Matrix:')
            print(analysis_operations.BinProbabilitiesFromRates(mean_rate_matrix)) 
        except:

            print('Singular rate matrix')
        print('Actual Bin probabilities:')
        print(analysis_operations.meanBinProbabilities(iterations, iteration_counter -10, iteration_counter))         

        
    # log iteration (Rainer)
    logger.log_iteration(iteration)
    
#logger.close()
print('hdWE sucessfully completed.                                           ')
