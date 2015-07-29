#!/usr/bin/python3
from __future__ import print_function
import numpy

from   lib.bin import Bin
import lib.constants as constants

class Iteration(object):
    """
    Defines one iteration of the hdWE algorithm
    """
    def __init__(self, iteration_id, boundaries, outer_region_boundaries):
        """ 
        @param ref_coords the path to reference coordinates defining the bin
        @param trajectories single or list of trajectories to 
               init the bin
        """
        # points to the reference structure of this bin
        self.iteration_id            = iteration_id    # int
        self.boundaries              = boundaries 
        self.outer_region_boundaries = outer_region_boundaries
        # the array of segments
        self.bins                    = []
        self.probability_flow        = 0.0                      

    def generateBin(self, reference_iteration_id, 
                    reference_bin_id, reference_segment_id,
                    target_number_of_segments, coordinate_ids,
                    outer_region):
        """
        Initialize a new instance of class Bin and append to bins
        @return  bin_id returns the id of the created bin
        """
        __bin = Bin(iteration_id               = self.getId(), 
                    bin_id                     = len(self.bins),  
                    reference_iteration_id     = reference_iteration_id, 
                    reference_bin_id           = reference_bin_id, 
                    reference_segment_id       = reference_segment_id, 
                    target_number_of_segments  = target_number_of_segments, 
                    coordinate_ids             = coordinate_ids,
                    outer_region               = outer_region)
        return self.__addBin(__bin)

    def __addBin(self, _bin):
        """
        Add the specified segments to this bin        
        @param
        """
        self.bins.append(_bin)
        return _bin.getId()
    
    def getId(self):
        return self.iteration_id
        
    def getBoundaries(self):
        return self.boundaries
    
    def getOuterRegionBoundaries(self):
        return self.outer_region_boundaries
        
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
    
    def getNumberOfActiveBins(self):
        """
        Number of bins
        """
        n = 0
        for this_bin in self.bins:
            if this_bin.outer_region == False:
                n += 1
        return n

    
    def getNumberOfSegments(self):
        """
        Returns the current number of segments
        """
        __number_of_segments = 0
        for __bin in self.bins:
            __number_of_segments += __bin.getNumberOfSegments()
        return __number_of_segments
    
    def getNumberOfActiveSegments(self):
        """
        Returns the current number of segments
        """
        __number_of_segments = 0
        for __bin in self.bins:
            if __bin.outer_region == False:
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
        Returns the target number of segments
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
        N = self.getNumberOfBins()
        flux_matrix = numpy.zeros((N,N), float)
        # Sum all probability that enters a bin from the parent bins
        for this_bin in self:
            for this_segment in this_bin.initial_segments:
                flux_matrix[this_segment.getParentBinId(), this_bin.getId()] += \
                                                this_segment.getProbability()
        return flux_matrix
        
    def RateMatrix(self):
        """
        Calculates the Rate Matrix for the given iteration.
        rate from bin i to bin j are stored in value [i][j]
        """
        # Initialize array for rate matrix
        N = self.getNumberOfBins()
        rate_matrix = numpy.zeros((N,N), float)
        # Sum all probability that enters a bin from the parent bins
        for this_bin in self:
            for this_segment in this_bin.initial_segments:
                rate_matrix[this_segment.getParentBinId(), this_segment.getBinId()] += \
                    sum(this_segment.getProbability())
        
        # Normalize all outrates with respect to a bin:
        for i in range(0,len(rate_matrix)):
            previous_iteration_parent_probability = 0.0
            # Sum up the parent bin's probabilities
            for j in range(0,len(rate_matrix)):
                previous_iteration_parent_probability += rate_matrix[i,j]
            if previous_iteration_parent_probability > constants.num_boundary:
                for j in range(0,len(rate_matrix)):
                    rate_matrix[i,j] /= previous_iteration_parent_probability
        return rate_matrix
    
    def TransitionMatrix(self):
        """
        Return the transition matrix for this iteration.
        Element [i,j] contains the number of segments migrating from
        bin i to bin j during this iteration. 
        """
        N = self.getNumberOfBins()
        transition_matrix = numpy.zeros((N,N), int)
        for this_bin in self:
            for this_segment in this_bin.initial_segments:
                transition_matrix[this_segment.getParentBinId(), this_segment.getBinId()] += 1
        return transition_matrix
    
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
    
    def getStartStateBins(self):
        """
        @return the list of bins in end state 
        """
        start_bins = []
        for this_bin in self.bins:
            if this_bin.isStartStateBin():
                start_bins.append(this_bin)
        return start_bins
    
    def getEndStateBins(self):
        """
        @return the list of bins in start state 
        """
        end_bins = []
        for this_bin in self.bins:
            if this_bin.isEndStateBin():
                end_bins.append(this_bin)
        return end_bins
    
    def setProbabilityFlow(self, probability_flow):
        self.probability_flow = probability_flow

    def getProbabilityFlow(self):
        return self.probability_flow
    
    def resetOuterRegion(self):
        """
        Remove the segments (not the initial_segments) from outer-region bins 
        and add the corresponding probability equally to the segments (no the initial_segments) of
        the parent bin
        """
        for this_bin in self:
            if this_bin.outer_region == True:
                for this_initial_segment in this_bin.initial_segments:
                    parent_bin_id = this_initial_segment.getParentBinId()
                    probability   = this_initial_segment.getProbability()
                    probability   = 1.0 * probability / self.bins[parent_bin_id].getNumberOfSegments()
                    for this_parent_segment in self.bins[parent_bin_id]:
                        #this_parent_segment.addProbability(probability)
                        # For some reason addProbability function does not work correctly here
                        this_parent_segment.probability = this_parent_segment.getProbability() + probability
                this_bin.deleteAllSegments()

        return
        
    def getOuterRegionFlag(self, coordinate_ids):
        """
        determines whether the bin lies in the outer region from information in the conf.file
        """

        for dimension in range(len(self.outer_region_boundaries)):
            if coordinate_ids[dimension] in self.outer_region_boundaries[dimension]:
                return True
   
        return False
            
            
        