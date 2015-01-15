class Segment(object):
    """
    Defines a trajectory element of a bin
    """
    def __init__(self, probability, parent_bin_id, parent_segment_id,
                 iteration_id, bin_id, segment_id):
        self.probability         = probability          # float
        self.parent_iteration_id = iteration_id - 1     # int
        self.parent_bin_id       = parent_bin_id        # int
        self.parent_segment_id   = parent_segment_id    # int
        self.bin_id              = bin_id               # int
        self.segment_id          = segment_id           # int
        self.iteration_id        = iteration_id         # int

    def getNameString(self):
        """Returns the indices in a string following the scheme iteration_bin_segment
        """
        return "{iteration:05d}-{_bin:05d}-{segment:05d}"\
        .format(iteration=self.iteration_id,
                _bin=self.bin_id,
                segment=self.segment_id)
    
    def getParentNameString(self):
        """Returns the indices in a string following the scheme iteration_bin_segment
        """
        return "{iteration:05d}-{_bin:05d}-{segment:05d}"\
                .format(iteration=self.parent_iteration_id,
                        _bin=self.parent_bin_id,
                        segment=self.parent_segment_id)
        
    def getProbability(self):
        return self.probability

    def addProbability(self, probability):
        self.probability += probability

    def subProbability(self, probability):
        self.probability -= probability

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
    
