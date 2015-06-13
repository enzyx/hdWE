import numpy as np

def getCoordinateIds(coordinates, bin_boundaries):
    """
    transforms coordinate values to coordinate bin ids
    segments exactly at the boundary are assigned to the bin with the smaller coordinate
    bin boundaries [1,2,3,4] lead to coordinate ids 0,1,2,3,4 where 0 is [0,1], 1 is ]1,2], 2is ]2,3], 3 is ]3,4] and 4 is ]4,inf
    """
    segment_coordinate_ids = []
    for dimension,coordinate in enumerate(coordinates):
        segment_coordinate_ids.append(-1)
        for dimension_bin_index, boundary in enumerate(bin_boundaries[dimension]):
            if coordinate <= boundary:
                segment_coordinate_ids[-1] = dimension_bin_index
                break
        if segment_coordinate_ids[-1] == -1:
            segment_coordinate_ids[-1] = len(bin_boundaries[dimension])
    
    return np.array(segment_coordinate_ids)