class Segment()
    """
    Defines a trajectory element of a bin
    """
    probability         = float(0.0)
    parent_iteration_id = int(0)
    parent_bin_id       = int(0)
    parent_segment_id   = int(0)
    iteration_id        = int(0)
    bin_id              = int(0)
    segment_id          = int(0)
    
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
        pass
    
    def getParentNameString(self):
        pass
        
    def getProbability(self):
        return self.probability
