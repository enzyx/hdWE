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
        """Returns the indices in a string following the scheme iteration_bin_segment
        """
        name_string= str(self.iteration_id).zfill(5) + '_' +str(self.bin_id).zfill(5) + '_' + str(self.segment_id).zfill(5) 
        return name_string
    
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
