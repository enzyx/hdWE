#!/usr/bin/python3
from __future__ import print_function
from bin import Bin
import numpy
import constants

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
        self.iteration_id = iteration_id    # int
        # the array of segments
        self.bins = []                      

    def generateBin(self, reference_iteration_id, 
                    reference_bin_id, reference_segment_id,
                    target_number_of_segments, outrates_converged = False):
        """
        Initialize a new instance of class Bin and append to bins
        @return  bin_id returns the id of the created bin
        """
        __bin = Bin(self.getId(), len(self.bins), reference_iteration_id, 
                    reference_bin_id, reference_segment_id,
                    target_number_of_segments, outrates_converged)
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
        
    def getNameString(self):
        """Returns iteration index as a string
        """
        name_string = str(self.iteration_id).zfill(5)
        return name_string

    def getProbability(self):
        """
        @returns The cumulative probability of all bins
        """
        probability = 0.0
        for _bin in self.bins:
            probability += _bin.getProbability()
        return probability
        
    def checkProbability(self):
        """
        @returns boolean: True if probability sums to 1
        """
        if self.getProbability() == 1.0:
            return True
        return False

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

    def getNumberOfPropagatedSegments(self):
        """
        Returns the current number of propagated(not converged) segments.
        """
        __number_of_segments = 0
        for __bin in self.bins:
            __number_of_segments += __bin.getNumberOfPropagatedSegments()
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

    def next(self):
        """
        Returns the next element of the array self.bins for pyton 2.6
        """
        self._iter_index += 1
        if self._iter_index >= len(self.bins):
            raise StopIteration
        else:
            return self.bins[self._iter_index]
            
    def __next__(self):
        """
        Returns the next element of the array self.bins
        """
        self._iter_index += 1
        if self._iter_index >= len(self.bins):
            raise StopIteration
        else:
            return self.bins[self._iter_index]
    
    def getSegments(self):
        """
        @return list of all segments in this iteration
        """
        segments = []
        for _bin in self.bins:
            for segment in _bin:
                segments.append(segment)
        return segments
            
    def FluxMatrix(self):
        """
        Calculates the Flux Matrix for the given iteration.
        """
        # Initialize array for flux matrix
        flux_matrix = numpy.zeros([self.getNumberOfBins(), self.getNumberOfBins()], float)
        # Sum all probability that enters a bin from the parent bins
        for bin_l in self.bins:
            for segment_l in bin_l.segments:
                flux_matrix[segment_l.getParentBinId(), bin_l.getId()] += \
                segment_l.getProbability()
        return flux_matrix
        
    def RateMatrix(self):
        """
        Calculates the Rate Matrix for the given iteration.
        rate from bin i to bin j are stored in value [i][j]
        """
        # Initialize array for rate matrix
        rate_matrix = numpy.zeros([self.getNumberOfBins(), self.getNumberOfBins()], float)
        # Sum all probability that enters a bin from the parent bins
        for bin_l in self.bins:
            for segment_l in bin_l.segments:
                if self.bins[segment_l.getParentBinId()].getProbability() > constants.num_boundary:
                    rate_matrix[segment_l.getParentBinId(), bin_l.getId()] += \
                    segment_l.getProbability()
        # Normalize all outrates with respect to a bin:
        for i in range(0,len(rate_matrix)):
            previous_iteration_parent_probability = 0.0            
            for j in range(0,len(rate_matrix)):
                previous_iteration_parent_probability += rate_matrix[i,j]
            if previous_iteration_parent_probability > constants.num_boundary:
                for j in range(0,len(rate_matrix)):
                    rate_matrix[i,j] /= previous_iteration_parent_probability
        return rate_matrix
    
    def getBinProbabilities(self):
        """
        Return a vector with the bin probablities
        """
        bin_probs = []
        for _bin in self.bins:
            bin_probs.append(_bin.getProbability())
        return bin_probs
    
    def getMaxBinProbability(self):
        """
        Returns the probablity of the bin with largest probability.
        """
        max_bin_probability = 0.0
        for bin_l in self.bins:
            buf_prob = bin_l.getProbability()
            if buf_prob > max_bin_probability:
                max_bin_probability = buf_prob

        return max_bin_probability        
