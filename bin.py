#!/usr/bin/python3

class Bin(object):
    """
    The bin class contains an array of segments (trajectories) and has a
    overall probability which is the sum of all segments probabilities
    """
    # points to the reference structure of this bin
    iteration_id              = int(0)
    bin_id                    = int(0)
    reference_iteration_id    = int(0)
    reference_bin_id          = int(0)
    reference_segment_id      = int(0)
    
    # the array of trajectories
    segments                  = []
    # How many segments we want in this bin:
    target_number_of_segments = int(0)
    
    def __init__(self, iteration_id, bin_id, reference_iteration_id, 
                 reference_bin_id, reference_segment_id, 
                 target_number_of_segments):
        """ 
        @param ref_coords the path to reference coordinates defining the bin
        @param trajectories single or list of trajectories to 
               init the bin
        """
        self.iteration_id              = iteration_id
        self.bin_id                    = bin_id
        self.reference_iteration_id    = reference_iteration_id
        self.reference_bin_id          = reference_bin_id
        self.reference_segment_id      = reference_segment_id
        self.target_number_of_segments = target_number_of_segments
    
    def generateSegment(self, probability, parent_bin_id, parent_segment_id):
        """
        @return segment_id of the created segment
        """
        __segment = Segment(probability=probability,
                              parent_bin_id=parent_bin_id,
                              parent_segment_id=parent_segment_id,
                              iteration_id=self.getIterationId(),
                              bin_id=self.getId(),
                              segment_id=len(self.segments))
        return self.__addSegment(__segment)
        
    def resampleSegments(self):
        """
        Split or Merge segments to generate the target number of segments
        """
        pass
    
    def __deleteSegments(self, segment_index):
        """
        @private should not be accessed from outside
        If the number of trajectories is smaller then N
        then we want to duplicated some trajectories 
        to fill the bin again
        """
        pass
    
    def __addSegment(self, segment):
        """
        @private should not be accessed from outside
        Add the specified segments to this bin        
        @param segment to add to segments
        @return the segment id of the added segment
        """
        self.segments.append(segment)
        return len(self.segments)

    def getStringName(self):
        pass

    def getId(self):
        """
        @return The id of this bin
        """
        return self.bin_id

    def getIterationId(self):
        """
        @return The iteration id this bin belongs to
        """
        return self.iteration_id
        
    def getReferenceIterationId(self):
        """
        @return The reference iteration id which points to the reference
                coordinates of this bin
        """
        return self.reference_iteration_id
    
    def getReferenceBinId(self):
        """
        @return The reference bin id which points to the reference
                coordinates of this bin
        """
        return self.reference_bin_id
    
    def getReferenceSegmentId(self):
        """
        @return The reference segment id which points to the reference
                coordinates of this bin
        """
        return self.reference_segment_id

    def getProbability(self):
        """
        Returns the cumulative probability of all bin trajectories
        """
        probability = 0.0
        for segment in self.segments:
            probability += segment.getProbability()
        return probability
    
    def getCurrentNumberOfSegments(self):
        """
        @return Current number of segments
        """
        return len(self.segments)
    
    def getTargetNumberOfSegments(self):
        """
        @return Target number of segments
        """
        return self.target_number_of_segments
    
    def __iter__(self):
        """
        Defines the class as iterable.
        """
        self._iter_index = -1
        return self

    def next(self):
        """
        Returns the next element of the array self.segments
        """
        self._iter_index += 1
        if self._iter_index >= len(self.segments):
            raise StopIteration
        else:
            return self.segments[self._iter_index]
