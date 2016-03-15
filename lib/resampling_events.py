
class Merge(object):
    """
    A class containing all the information of a merge event
    
    surviving_segment: The index of the surviving segment
    deleted_segments: List of segment_ids which are 
                      deleted probability
    """
    def __init__(self, surviving_segment, deleted_segments):
        self.surviving_segment = surviving_segment
        self.deleted_segments  = deleted_segments
    
    def getType(self):
        return type(self).__name__
    
class Split(object):
    """
    A class containing all the information of a split event
    
    parent_segment: segment_id of the split segment
    m: Number of segment which result from splitting
    """
    def __init__(self, parent_segment, m):
        self.parent_segment  = parent_segment
        self.m               = m
    
    def getType(self):
        return type(self).__name__
        