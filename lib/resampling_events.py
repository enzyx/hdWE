
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
        