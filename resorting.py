import numpy

def copyBinStructureToLastIteration(iterations, SEGMENTS_PER_BIN):
    """
    generates new empty bins in the newest iteration for all bins of previous iteration 
    """
    current_iteration = iterations[-1]
    parent_iteration  = iterations[-2]
    
    for parent_bin in parent_iteration:
        current_iteration.generateBin(reference_iteration_id    = parent_bin.getReferenceIterationId(),
                                      reference_bin_id          = parent_bin.getReferenceBinId(),
                                      reference_segment_id      = parent_bin.getReferenceSegmentId(),
                                      target_number_of_segments = SEGMENTS_PER_BIN,
                                      rama_id                   = parent_bin.getRamaId(),
                                      outrates_converged        = parent_bin.isConverged())

def resort(iterations, md_module, COORDINATE_THRESHOLD, SEGMENTS_PER_BIN, MAX_NUMBER_OF_BINS):
    """
    Resorts final segments of parent iteration into bins and creates new bins if necessary
    """
    parent_iteration    = iterations[-2]
    current_iteration   = iterations[-1]
    rmsd_matrix         = md_module.calcRmsdSegmentsToBinsMatrix(parent_iteration)
    
    # Sanity check for RMSD matrix
#     if rmsd_matrix.min() == 0.0:
#         print("\x1b[31m Warning!\x1b[0m SegmentsToBins RMSD Matrix has entries 0.00000!"\
#               " This is extremely unlikely!")
    # Save the rmsd matrix there for checking...
    numpy.savetxt("{0}-log/{1}.dihe".format(md_module.jobname, current_iteration.getNameString()), 
                  rmsd_matrix,
                  fmt='%d')
    
    segment_id = -1
    for parent_bin in parent_iteration:
        for parent_segment in parent_bin:
            is_segment_handled = False
            segment_id += 1
            # 1. Check if parent_segment fits into a bin of
            #    the previous iteration
            for this_bin in current_iteration.bins:
                if rmsd_matrix[segment_id] == this_bin.getRamaId():
                    this_bin.generateSegment(probability         = parent_segment.getProbability(),
                                             parent_iteration_id = parent_segment.getIterationId(),
                                             parent_bin_id       = parent_segment.getBinId(),
                                             parent_segment_id   = parent_segment.getId())
                    is_segment_handled = True
                    break
                                   
            # 2. If it fits nowhere create new bin
            if not is_segment_handled and not current_iteration.getNumberOfBins() >= MAX_NUMBER_OF_BINS: 
                bin_id = current_iteration.generateBin(reference_iteration_id      = parent_segment.getIterationId(),
                                                       reference_bin_id            = parent_segment.getBinId(),
                                                       reference_segment_id        = parent_segment.getId(),
                                                       target_number_of_segments   = SEGMENTS_PER_BIN,
                                                       rama_id                     = rmsd_matrix[segment_id])
                current_iteration.bins[bin_id].generateSegment(probability         = parent_segment.getProbability(),
                                                               parent_iteration_id = parent_segment.getIterationId(),
                                                               parent_bin_id       = parent_segment.getBinId(),
                                                               parent_segment_id   = parent_segment.getId())
                is_segment_handled = True
            
            # 3. We were not allowed to create a new bin
            #    so stuff it into its parent bin
            if not is_segment_handled:
                current_iteration.bins[parent_segment.getBinId()].generateSegment(
                                             probability         = parent_segment.getProbability(),
                                             parent_iteration_id = parent_segment.getIterationId(),
                                             parent_bin_id       = parent_segment.getBinId(),
                                             parent_segment_id   = parent_segment.getId())
                is_segment_handled = True
            # 4. Sanity check
            if not is_segment_handled:
                raise("Somthin' smells fishy. Segment was not sorted into any bin.")
    return