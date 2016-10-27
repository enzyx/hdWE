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
def copyBinStructureToLastIteration(iterations):
    """
    generates new empty bins in the newest iteration for all bins of previous iteration 
    """
    current_iteration = iterations[-1]
    parent_iteration  = iterations[-2]
    
    for parent_bin in parent_iteration:
        current_iteration.generateBin(target_number_of_segments = parent_bin.getTargetNumberOfSegments(),
                                      coordinate_ids            = parent_bin.getCoordinateIds(),
                                      sample_region             = parent_bin.getSampleRegion()) 

def resort(iterations, md_module, INITIAL_TARGET_NUMBER_OF_SEGMENTS):
    """
    Resorts final segments of parent iteration into bins and creates new bins if necessary
    """
    parent_iteration    = iterations[-2]
    current_iteration   = iterations[-1]
    
    # Sanity check for RMSD matrix
    #     if rama_ids.min() == 0.0:
    #         print("\x1b[31m Warning!\x1b[0m SegmentsToBins RMSD Matrix has entries 0.00000!"\
    #               " This is extremely unlikely!")
    
    for parent_bin in parent_iteration:
        for parent_segment in parent_bin:
            is_segment_handled = False
            these_coordinate_ids = parent_segment.getCoordinateIds(current_iteration.getBoundaries())
                
            # 1. Check if parent_segment fits into a bin of
            #    the previous iteration
            if not is_segment_handled:
                for this_bin in current_iteration.bins:
                    if all(these_coordinate_ids == this_bin.getCoordinateIds()):
                        this_bin.generateSegment(probability         = parent_segment.getProbability(),
                                                 parent_iteration_id = parent_segment.getIterationId(),
                                                 parent_bin_id       = parent_segment.getBinId(),
                                                 parent_segment_id   = parent_segment.getId())
                        is_segment_handled = True
                        break
                                   
            # 2. If it fits nowhere create new bin
            if not is_segment_handled: 
                bin_id = current_iteration.generateBin(target_number_of_segments   = INITIAL_TARGET_NUMBER_OF_SEGMENTS,
                                                       coordinate_ids              = these_coordinate_ids,
                                                       sample_region               = current_iteration.isInSampleRegion(these_coordinate_ids))
                current_iteration.bins[bin_id].generateSegment(probability         = parent_segment.getProbability(),
                                                               parent_iteration_id = parent_segment.getIterationId(),
                                                               parent_bin_id       = parent_segment.getBinId(),
                                                               parent_segment_id   = parent_segment.getId())
                is_segment_handled = True
            
            # 3. Sanity check
            if not is_segment_handled:
                raise("Somthin' smells fishy. Segment was not sorted into any bin.")
            
