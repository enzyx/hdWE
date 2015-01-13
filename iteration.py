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

    def generateBin(self, reference_iteration_id, 
                    reference_bin_id, reference_segment_id, target_number_of_segments)
        """
        Initialize a new instance of class Bin and append to bins
        @return  bin_id returns the id of the created bin
        """
        tmp_bin = Bin(self.getId(), len(self.bins), reference_iteration_id, 
                      reference_bin_id, reference_segment_id, target_number_of_segments)
        return self._addBin(tmp_bin)

    def __addBin(self, _bin):
        """
        Add the specified segments to this bin        
        @param
        """
        self.bins.append(_bin)
        return len(self.bins)
    
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
        return len(segments)
    
    def getNumberOfSegments(self):
        pass
