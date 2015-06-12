import numpy
import bin_classifier      

def copyBinStructureToLastIteration(iterations):
    """
    generates new empty bins in the newest iteration for all bins of previous iteration 
    """
    current_iteration = iterations[-1]
    parent_iteration  = iterations[-2]
    
    for parent_bin in parent_iteration:
        current_iteration.generateBin(reference_iteration_id    = parent_bin.getReferenceIterationId(),
                                      reference_bin_id          = parent_bin.getReferenceBinId(),
                                      reference_segment_id      = parent_bin.getReferenceSegmentId(),
                                      target_number_of_segments = parent_bin.getTargetNumberOfSegments(),
                                      coordinate_ids            = parent_bin.getCoordinateIds(),
                                      start_states              = parent_bin.isStartStateBin(),
                                      end_states                = parent_bin.isEndStateBin())

def resort(iterations, md_module, INITIAL_TARGET_NUMBER_OF_SEGMENTS, START_STATES, END_STATES, STEADY_STATE):
    """
    Resorts final segments of parent iteration into bins and creates new bins if necessary
    """
    parent_iteration    = iterations[-2]
    current_iteration   = iterations[-1]
    coordinate_ids      = md_module.calcCoordinateIds(parent_iteration)
    
    # Sanity check for RMSD matrix
    #     if rama_ids.min() == 0.0:
    #         print("\x1b[31m Warning!\x1b[0m SegmentsToBins RMSD Matrix has entries 0.00000!"\
    #               " This is extremely unlikely!")
    
    segment_id = -1
    probability_flow = 0.0
    for parent_bin in parent_iteration:
        for parent_segment in parent_bin:
            is_segment_handled = False
            segment_id += 1
            
            # 1. Check if parent_segment is in an end state bin
            if STEADY_STATE:
                for end_state_bin in parent_iteration.getEndStateBins():
                    if coordinate_ids[segment_id] == end_state_bin.getCoordinateIds():
                        probability_flow += parent_segment.getProbability()
                        is_segment_handled = True
                        break
                
            # 2. Check if parent_segment fits into a bin of
            #    the previous iteration
            if not is_segment_handled:
                for this_bin in current_iteration.bins:
                    if coordinate_ids[segment_id] == this_bin.getCoordinateIds():
                        this_bin.generateSegment(probability         = parent_segment.getProbability(),
                                                 parent_iteration_id = parent_segment.getIterationId(),
                                                 parent_bin_id       = parent_segment.getBinId(),
                                                 parent_segment_id   = parent_segment.getId())
                        is_segment_handled = True
                        break
                                   
            # 3. If it fits nowhere create new bin
            if not is_segment_handled: 
                bin_id = current_iteration.generateBin(reference_iteration_id      = parent_segment.getIterationId(),
                                                       reference_bin_id            = parent_segment.getBinId(),
                                                       reference_segment_id        = parent_segment.getId(),
                                                       target_number_of_segments   = INITIAL_TARGET_NUMBER_OF_SEGMENTS,
                                                       coordinate_ids              = coordinate_ids[segment_id],
                                                       start_states                = START_STATES,
                                                       end_states                  = END_STATES)
                current_iteration.bins[bin_id].generateSegment(probability         = parent_segment.getProbability(),
                                                               parent_iteration_id = parent_segment.getIterationId(),
                                                               parent_bin_id       = parent_segment.getBinId(),
                                                               parent_segment_id   = parent_segment.getId())
                is_segment_handled = True
            
            # 4. Sanity check
            if not is_segment_handled:
                raise("Somthin' smells fishy. Segment was not sorted into any bin.")
            
    # reassign flown probility
    print("Probability flow: {0:f}".format(probability_flow))                        
    current_iteration.setProbabilityFlow(probability_flow)         
    

def assignRecycledProbability(iteration):
        
    if iteration.getProbabilityFlow() > 0.0:
        # Check if segments exist in start state bins
        start_segments = 0
        # It is assumed that the starting structure lies in the starting state 
        # then the corresponding bin id is 0
        for start_bin in iteration.getStartStateBins():
            start_segments += start_bin.getNumberOfSegments()
        
        # No segments exist in start state bins, respawn a segment
        if start_segments == 0:
            iteration.bins[0].generateSegment(probability         = iteration.getProbabilityFlow(),
                                              parent_iteration_id = 0,
                                              parent_bin_id       = 0,
                                              parent_segment_id   = 0)           
        
        # Segments exist: Scale all segment probabilities accordingly
        else:
            # get normalization factor
            start_state_probability = 0.0
            for start_bin in iteration.getStartStateBins():
                start_state_probability += start_bin.getProbability()
    
            scale_factor = 1.0 + iteration.getProbabilityFlow() / start_state_probability
            
            # scale probabilities of START_STATE segments
            for start_bin in iteration.getStartStateBins():
                for this_segment in start_bin:
                        this_segment.multProbability(scale_factor)   
                        