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
    
    def newBin()
        """
        Initialize a new instance of class Bin and append to bins
        """
        #tmp_bin = Bin(a,b,c,d)
        #self.addBin(tmp_bin)

    def addBin(self, _bin):
        """
        Add the specified segments to this bin        
        @param
        """
        self.bins.append(_bin)


