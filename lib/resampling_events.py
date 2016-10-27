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
class Merge(object):
    """
    A class containing all the information of a merge event
    
    surviving_segment: The index of the surviving segment
    deleted_segments: List of segment_ids which are 
                      deleted probability
    """
    def __init__(self, surviving_segment_id, deleted_segments_ids):
        self.surviving_segment_id  = surviving_segment_id
        self.deleted_segments_ids  = deleted_segments_ids
    
    def getType(self):
        return type(self).__name__
    
class Split(object):
    """
    A class containing all the information of a split event
    
    parent_segment: segment_id of the split segment
    m: Number of segment which result from splitting
    """
    def __init__(self, parent_segment_id, m):
        self.parent_segment_id  = parent_segment_id
        self.m                  = m
    
    def getType(self):
        return type(self).__name__
        
