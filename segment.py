class Segment(object):
    """
    Defines a trajectory element of a bin
    """
    def __init__(self, probability, parent_bin_id, parent_segment_id,
                 iteration_id, bin_id, segment_id):
        self.probability         = probability
        self.parent_bin_id       = parent_bin_id
        self.parent_segment_id   = parent_segment_id
        self.bin_id              = bin_id
        self.segment_id          = segment_id
        self.iteration_id        = iteration_id
        self.parent_iteration_id = iteration_id - 1

    def getNameString(self):
        return "{iteration:05d}-{_bin:05d}-{segment:05d}"\
        .format(iteration=self.iteration_id,
                _bin=self.bin_id,
                segment=self.segment_id)
    
    def getParentNameString(self):
        return "{iteration:05d}-{_bin:05d}-{segment:05d}"\
                .format(iteration=self.parent_iteration_id,
                        _bin=self.parent_bin_id,
                        segment=self.parent_segment_id)
        
    def getProbability(self):
        return self.probability

    def getId(self):
        return self.segment_id
    
    def getBinId(self):
        return self.bin_id

    def getIterationId(self):
        return self.iteration_id

    def getParentIterationId(self):
        return self.parent_iteration_id   
                
    def getParentBinId(self):
        return self.parent_bin_id
        
    def getParentSegmentId(self):
        return self.parent_segment_id
    
