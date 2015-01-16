#!/usr/bin/python3

import sys
import os
import arg_parser as args
import numpy
import initiate
from iteration import Iteration
from bin import Bin
from segment import Segment
from logger import Logger

# The global list of arrays
iterations = [] 

# Initialize the logger
logger = Logger(args.work_dir+args.logfile)
# Read previous log
if if os.stat(logfile).st_size != 0:
    iterations = logger.load_iterations()

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
                
    if debug: 
        print("The overall probabiliy is at {0:05f}".format(iteration.getProbability()))
    
    # Run MD
    md_module.RunMDs(iteration)
    
    iterations.append(iteration) 
    
    # log iteration (Rainer)
    logger.log_iteration(iteration)
    
logger.close()
