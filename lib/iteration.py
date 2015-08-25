#!/usr/bin/python3
from __future__ import print_function
import numpy

from   lib.bin import Bin
import lib.constants as constants

class Iteration(object):
    """
    Defines one iteration of the hdWE algorithm
    """
    def __init__(self, iteration_id, boundaries, sample_region):
        """ 
        @param ref_coords the path to reference coordinates defining the bin
        @param trajectories single or list of trajectories to 
               init the bin
        """
        # points to the reference structure of this bin
        self.iteration_id            = iteration_id    # int
        self.boundaries              = boundaries 
        self.sample_region           = sample_region
        # the array of segments
        self.bins                    = []
        self.probability_flow        = 0.0                      

    def generateBin(self, target_number_of_segments, coordinate_ids,
                    sample_region):
        """
        Initialize a new instance of class Bin and append to bins
        @return  bin_id returns the id of the created bin
        """
        __bin = Bin(iteration_id               = self.getId(), 
                    bin_id                     = len(self.bins),  
                    target_number_of_segments  = target_number_of_segments, 
                    coordinate_ids             = coordinate_ids,
                    sample_region              = sample_region)
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
    
    def getSampleRegion(self):
        return self.sample_region
        
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
            if this_bin.sample_region == True:
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
            if __bin.sample_region == True:
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
        for this_bin in self.bins:
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
        for this_bin in self.bins:
            for this_segment in this_bin.initial_segments:
                # Check for numpy arrays and sum if required
                tot_segment_prob = 0.0
                # Some probability values are numpy.float64 type?
                if type(this_segment.getProbability()) != float \
                 and type(this_segment.getProbability()) != numpy.float64:
                    tot_segment_prob = sum(this_segment.getProbability())
                else:
                    tot_segment_prob = this_segment.getProbability()
                rate_matrix[this_segment.getParentBinId(), this_segment.getBinId()] += \
                    tot_segment_prob
        
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
        for this_bin in self.bins:
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
        for this_bin in self.bins:
            if this_bin.sample_region == False:
                for this_initial_segment in this_bin.initial_segments:
                    parent_bin_id = this_initial_segment.getParentBinId()
                    probability   = this_initial_segment.getProbability()
                    probability   = 1.0 * probability / self.bins[parent_bin_id].getNumberOfSegments()
                    if self.bins[parent_bin_id].getNumberOfSegments() > 0:
                        for this_parent_segment in self.bins[parent_bin_id]:
                            # this_parent_segment.addProbability(probability)
                            # For some reason addProbability function does not work correctly here
                            this_parent_segment.probability = this_parent_segment.getProbability() + probability
                    else:
                        print('Warning: Probability flow in to outer region could not be reset to parent bin {}, \
                        because parent bin is empty. Probability is deleted'.format(parent_bin = parent_bin_id))
                this_bin.deleteAllSegments()

        return
        
    def isInSampleRegion(self, coordinate_ids):
        """
        determines whether the bin lies in the outer region from information in the conf.file
        the region given in the conf.file masks the inner region bins,
        boundary 0 including
        boundary 1 excluding
        """
        dimension_flag = []
        for dimension in range(len(self.sample_region)):
            if coordinate_ids[dimension] == 0:
                bin_boundary_0 = - 1e99
                bin_boundary_1 = self.boundaries[dimension][0]
            elif coordinate_ids[dimension] == len(self.boundaries[dimension]):
                bin_boundary_0 = self.boundaries[dimension][-1]            
                bin_boundary_1 =  1e99
            else:
                bin_boundary_0 = self.boundaries[dimension][coordinate_ids[dimension] - 1 ]
                bin_boundary_1 = self.boundaries[dimension][coordinate_ids[dimension]]
                
            if self.sample_region[dimension][0] <= bin_boundary_1 and \
               self.sample_region[dimension][1] > bin_boundary_0:
                dimension_flag.append(True)
            else:
                dimension_flag.append(False)

        for dimension in dimension_flag:
            if dimension == False:
                return False
   
        return True
    
    def getNumberOfEmptyBins(self):
        """
        returns the number of current empty active bins
        """
        empty_bins = 0
        for bin_loop_tmp in self.bins:
            if bin_loop_tmp.getNumberOfSegments() == 0 and bin_loop_tmp.getSampleRegion() == True:
                empty_bins += 1
        return empty_bins