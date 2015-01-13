class Segment(object):
    """
    Defines a trajectory element of a bin
    """
    def __init__(self, probability, parent_bin_id, parent_segment_id,
                 iteration_id, bin_id, segment_id):
        self.probability         = probability
        self.parent_iteration_id = iteration_id - 1
        self.parent_bin_id       = parent_bin_id
        self.parent_segment_id   = parent_segment_id
        self.bin_id              = bin_id
        self.segment_id          = segment_id
        self.iteration_id        = iteration_id

    def getNameString(self):
        pass
    
    def getParentNameString(self):
        pass
        
    def getProbability(self):
        return self.probability

    def getId(self):
        return self.segment_id
    
    def getBinId(self):
        return self.bin_id

    def getIterationId(self):
        return self.iteration_id
    
    def getParentIterationId(self):
        return self.parent_iteration
    
    def getParentSegmentId(self):
        return self.parent_segment_id
    
    def getParentBinId(self):
        return self.parent_bin_id
