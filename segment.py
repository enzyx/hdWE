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

    def __getNameString(self, iteration_id, bin_id, segment_id):
        """
        @return Formated segment file name string
        """
        return "{iteration:05d}_{_bin:05d}_{segment:05d}"\
                .format(iteration=iteration_id,
                       _bin=bin_id,
                       segment=segment_id)
    
    def setSegmentId(self, segment_id):
        """
        required for resampling function to reorder the segment ids. 
        Should not be used from outside Bin class
        """
        self.segment_id = segment_id
        
    def setParentIterationId(self, parent_iteration_id):
        """
        required for reweighting function to create new segments in empy bins. 
        Should not be used from Bin class
        """ 
        self.parent_iteration_id = parent_iteration_id

    def getNameString(self):
        """
        @return the indices in a string following the scheme iteration_bin_segment
        """
        return self.__getNameString(iteration_id=self.iteration_id,
                                    bin_id=self.bin_id,
                                    segment_id=self.segment_id)
    
    def getParentNameString(self):
        """
        @return the indices in a string following the scheme iteration_bin_segment
        """
        return self.__getNameString(iteration_id=self.parent_iteration_id,
                                    bin_id=self.parent_bin_id,
                                    segment_id=self.parent_segment_id)
        
    def getProbability(self):
        return self.probability
        
    def setProbability(self, probability):
        self.probability = probability

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
        
    def __eq__(self, other_segment): 
        return self.__dict__ == other_segment.__dict__
    
    def isParent(self, segment):
        if segment.getParentIterationId() == self.getIterationId() and\
           segment.getParentBinId() == self.getBinId() and\
           segment.getParentSegmentId() == self.getId():
            return True
        return False

    
