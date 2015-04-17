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
                                      outrates_converged        = parent_bin.isConverged())

def resort(iterations, md_module, COORDINATE_THRESHOLD, SEGMENTS_PER_BIN):
    """
    Resorts final segments of parent iteration into bins and creates new bins if necessary
    """
    parent_iteration    = iterations[-2]
    current_iteration   = iterations[-1]
    #TODO: check rmsd_matrix calculation (i<->j?)
    rmsd_matrix         = md_module.calcRmsdSegmentsToBinsMatrix(parent_iteration)

    new_bins   = []
    for parent_bin in parent_iteration:
        for parent_segment in parent_bin:
            is_segment_handled = False
     
            # 1. Check if parent_segment fits into a bin of
            #    the previous iteration
            for bin_id, coordinate in enumerate(rmsd_matrix[parent_segment.getId(),:]):
                if coordinate <= COORDINATE_THRESHOLD and not is_segment_handled:
                    current_iteration.bins[bin_id].generateSegment(probability       = parent_segment.getProbability(),
                                                                   parent_bin_id     = parent_segment.getBinId(),
                                                                   parent_segment_id = parent_segment.getId())
                    is_segment_handled = True
                    break

            # 2. Check if parent_segment fits into one of the
            #    bins newly created for this iteration
            if not is_segment_handled and len(new_bins) > 0:
                new_rmsds = md_module.calcRmsdToBins(parent_segment, new_bins)
                for newbin_id, coordinate in enumerate(new_rmsds):
                    if coordinate <= COORDINATE_THRESHOLD:
                        current_iteration.bins[new_bins[newbin_id].getId()].\
                            generateSegment(probability       = parent_segment.getProbability(),
                                            parent_bin_id     = parent_segment.getBinId(),
                                            parent_segment_id = parent_segment.getId())
                        is_segment_handled = True                        
                        break
                                   
            # 3. If it fits nowhere create new bin            
            if not is_segment_handled:
                bin_id = current_iteration.generateBin(reference_iteration_id    = parent_segment.getIterationId(),
                                                       reference_bin_id          = parent_segment.getBinId(),
                                                       reference_segment_id      = parent_segment.getId(),
                                                       target_number_of_segments = SEGMENTS_PER_BIN)
                current_iteration.bins[bin_id].generateSegment(probability       = parent_segment.getProbability(),
                                                               parent_bin_id     = parent_segment.getBinId(),
                                                               parent_segment_id = parent_segment.getId())
                new_bins.append(current_iteration.bins[bin_id])
                is_segment_handled = True
                
            # 4. Sanity check
            if not is_segment_handled:
                raise("Somthin' smells fishy. Segment was not sorted into any bin.")
    return