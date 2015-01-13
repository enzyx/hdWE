#!/usr/bin/python3
from bin import Bin
from segment import Segment

class Iteration(object):
    """
    Defines one iteration of the hdWE algorithm
    """
    def __init__(self, iteration_id):
        """ 
        @param ref_coords the path to reference coordinates defining the bin
        @param trajectories single or list of trajectories to 
               init the bin
        """
        # points to the reference structure of this bin
        self.iteration_id = iteration_id
        # the array of segments
        self.bins = []

    def generateBin(self, reference_iteration_id, 
                    reference_bin_id, reference_segment_id,
                    target_number_of_segments):
        """
        Initialize a new instance of class Bin and append to bins
        @return  bin_id returns the id of the created bin
        """
        __bin = Bin(self.getId(), len(self.bins), reference_iteration_id, 
                    reference_bin_id, reference_segment_id,
                    target_number_of_segments)
        return self.__addBin(__bin)

    def __addBin(self, _bin):
        """
        Add the specified segments to this bin        
        @param
        """
        self.bins.append(_bin)
        return len(self.bins)-1
    
    def getId(self):
        return self.iteration_id

    def getProbability(self):
        """
        @returns The cumulative probability of all bins
        """
        probability = 0.0
        for _bin in self.bins:
            probability += _bin.getProbability()
        return probability

    def getNumberOfBins(self):
        """
        Number of bins
        """
        return len(self.bins)
    
    def getNumberOfSegments(self):
        """
        Returns the current number of segments
        """
        __number_of_segments = 0
        for __bin in self.bins:
            __number_of_segments += __bin.getNumberOfSegments()
        return __number_of_segments

    def getTargetNumberOfSegments(self):
        """
        Returns the current number of segments
        """
        __number_of_segments = 0
        for __bin in self.bins:
            __number_of_segments += __bin.getTargetNumberOfSegments()
        return __number_of_segments

    def __iter__(self):
        """
        Defines the class as iterable.
        """
        self._iter_index = -1
        return self

    def __next__(self):
        """
        Returns the next element of the array self.bins
        """
        self._iter_index += 1
        if self._iter_index >= len(self.bins):
            raise StopIteration
        else:
            return self.bins[self._iter_index]
