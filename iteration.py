#!/usr/bin/python3

class Iteration(object):
    # points to the reference structure of this bin
    iteration_id = int(0)
    # the array of trajectories
    bins         = []

    def __init__(self, iteration_id):
        """ 
        @param ref_coords the path to reference coordinates defining the bin
        @param trajectories single or list of trajectories to 
               init the bin
        """
        self.iteration_id = iteration_id

    def getProbability(self):
        """
        Returns the cumulative probability of all bin trajectories
        """
        pass
    
    def getNumberOfBins(self):
        """
        Number of bins
        """
        return len(segments)
    
    def getNumberOfSegments(self):
        pass
    
    def generateBin(self, reference_iteration_id, 
                    reference_bin_id, reference_segment_id)
        """
        Initialize a new instance of class Bin and append to bins
        @return  bin_id returns the id of the created bin
        """
        tmp_bin = Bin(self.getId(), len(self.bins), reference_iteration_id, 
                      reference_bin_id, reference_segment_id)
        self._addBin(tmp_bin)
        return len(self.bins)

    def _addBin(self, _bin):
        """
        Add the specified segments to this bin        
        @param
        """
        self.bins.append(_bin)
    
    def getId(self):
        return self.iteration_id
