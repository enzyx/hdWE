#!/usr/bin/python3

class Bin(object):
    # points to the reference structure of this bin
    iteration_id            = int(0)
    bin_id                  = int(0)
    reference_iteration_id  = int(0)
    reference_bin_id        = int(0)
    reference_segment_id    = int(0)
    
    # the array of trajectories
    segments                  = []
    # How many segments we want in this bin:
    target_number_of_segments = int(0)
    
    def __init__(self, iteration_id, bin_id, reference_iteration_id, 
                 reference_bin_id, reference_segment_id):
        """ 
        @param ref_coords the path to reference coordinates defining the bin
        @param trajectories single or list of trajectories to 
               init the bin
        """
        self.iteration_id           = iteration_id
        self.bin_id                 = bin_id
        self.reference_iteration_id = reference_iteration_id
        self.reference_bin_id       = reference_bin_id
        self.reference_segment_id   = reference_segment_id
    
    def generateSegment(self, probability, parent_bin_id, parent_segment_id):
        """
        @return segment_id of the created segment
        """
        tmp_segment = Segment(probability=probability,
                              parent_bin_id=parent_bin_id, 
                              parent_segment_id=parent_segment_id,
                              iteration_id=self.getIterationId(),
                              bin_id=self.getId(),
                              segment_id=len(self.segments))
        self._addSegment(tmp_segment)
        return len(self.segments)
        
    def getProbability(self):
        """
        Returns the cumulative probability of all bin trajectories
        """
        pass
    
    def getNumberOfCurrentSegments(self):
        """
        Number of trajectories
        """
        return len(segments)
    
    def resampleSegments(self):
        """
        Split or Merge segments to generate the target number of segments
        """
        pass
    
    def _deleteSegments(self, segment_index):
        """
        If the number of trajectories is smaller then N
        then we want to duplicated some trajectories 
        to fill the bin again
        """
        pass
    
    def _addSegment(self, segment):
        """
        Add the specified segments to this bin        
        @param
        """
        self.segments.append(segments)

    def getStringName(self):
        pass

    def getId(self):
        return self.bin_id

    def getIterationId(self):
        return self.iteration_id
        
    def getReferenceIterationId(self):
        return self.reference_iteration_id
    
    def getReferenceBinId(self):
        return self.reference_bin_id
    
    def getReferenceSegmentId(self):
        return self.reference_segment_id
