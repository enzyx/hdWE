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
