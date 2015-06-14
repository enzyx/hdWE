#!/usr/bin/python3
from __future__ import print_function
from segment import Segment
import random as rnd
import copy

class Bin(object):
    """
    The bin class contains an array of segments (trajectories) and has a
    overall probability which is the sum of all segments probabilities
    """
    def __init__(self, iteration_id, bin_id, reference_iteration_id, 
                 reference_bin_id, reference_segment_id, 
                 target_number_of_segments, coordinate_ids, start_states, end_states):
        """
        @param ref_coords the path to reference coordinates defining the bin
        @param trajectories single or list of trajectories to 
               init the bin
        """
        # points to the reference structure of this bin
        self.iteration_id              = iteration_id               # int
        self.bin_id                    = bin_id                     # int
        self.reference_iteration_id    = reference_iteration_id     # int
        self.reference_bin_id          = reference_bin_id           # int
        self.reference_segment_id      = reference_segment_id       # int
        # How many segments we want in this bin
        self.target_number_of_segments = target_number_of_segments  # int
        # coordinate bin ids for considered binning dimensions
        self.coordinate_ids            = coordinate_ids
        self.is_start_bin              = self.__testBinState(start_states)
        self.is_end_bin                = self.__testBinState(end_states) 
        # The array of segments
        self.segments                  = []
        # In this array the segments are copied before resampling happens.
        # We need to store the old segments information
        # to be able to correctly recalculate the bin to bin rates 
        # after resampling.
        self.initial_segments          = []

    def generateSegment(self, probability, parent_iteration_id, parent_bin_id, parent_segment_id):
        """
        @return segment_id of the created segment
        """
        __segment = Segment(probability         = probability,
                            parent_iteration_id = parent_iteration_id,
                            parent_bin_id       = parent_bin_id,
                            parent_segment_id   = parent_segment_id,
                            iteration_id        = self.getIterationId(),
                            bin_id              = self.getId(),
                            segment_id          = len(self.segments))
        return self.__addSegment(__segment)

    def respawnSegmentFromReference(self, probability):
        __segment = self.generateSegment(
                             probability          = probability,
                             parent_iteration_id  = self.getReferenceIterationId(),
                             parent_bin_id        = self.getReferenceBinId(),
                             parent_segment_id    = self.getReferenceSegmentId())
        self.setConverged(False)
        return __segment
        
    def resampleSegments(self, merge_mode=0):
        """
        Split or Merge segments to generate the target number of segments
        """
        if len(self.segments) == 0:
            return 
        # Too many bins -> merge
        prob_tot = self.getProbability()
        if len(self.segments) > self.target_number_of_segments:
            #Merge according to segment probabilities            
            if merge_mode == 0:
                for c in range(len(self.segments) - self.target_number_of_segments):
                    # Get the extinction index
                    ext_index = 0
                    merge_index = 0
                    inv_weights = []
                    inv_weights_tot = 0.0
                    for segment in self.segments:
                        inv_weights_tot += prob_tot/segment.getProbability()
                        inv_weights.append(prob_tot/segment.getProbability())
                    extinction_probabilities = []
                    for inv_weight in inv_weights:
                        extinction_probabilities.append(inv_weight/inv_weights_tot)
                    rand = rnd.random()
                    cumulated_probability = 0.0
                    for index in range(len(extinction_probabilities)):
                        cumulated_probability += extinction_probabilities[index]
                        if cumulated_probability >= rand:
                            ext_index = index
                            break
                    # Get now the merge index 
                    inv_weights = []
                    inv_weights_tot = 0.0
                    for index in range(len(self.segments)):
                        if index == ext_index:
                            inv_weights.append(0.0)
                            continue
                        inv_weights_tot += prob_tot/self.segments[index].getProbability()
                        inv_weights.append(prob_tot/self.segments[index].getProbability())
                    extinction_probabilities = []
                    for inv_weight in inv_weights:
                        extinction_probabilities.append(inv_weight/inv_weights_tot)
                    rand = rnd.random()
                    cumulated_probability = 0.0
                    for index in range(len(extinction_probabilities)):
                        cumulated_probability += extinction_probabilities[index]
                        if cumulated_probability >= rand:
                            merge_index = index
                            break
                    # Pew, now we have to indices, a merge and a extinction index
                    # Do the merging stuff now.
                    shift_prob = self.segments[ext_index].getProbability()
                    self.segments[merge_index].addProbability(shift_prob)
                    del self.segments[ext_index]
                    # Reorder segment ids after deletion 
                    self.__fixSegmentIds()
                    #TODO: reset the indices 
                return
            #Random Merge ignoring segment probabilities
            elif merge_mode == 1:
                init_num_segments = len(self.segments)
                del_prob = 0.0
                for i in range(init_num_segments,self.target_number_of_segments, -1):
                    #print(str(i))
                    ext_index = rnd.randint(0, i - 1 )
                    #print(str(ext_index))
                    del_prob += self.segments[ext_index].getProbability()
                    del self.segments[ext_index]
                    self.__fixSegmentIds()

                del_prob = 1.0 * del_prob / self.getNumberOfSegments()                
                for i in range(0,self.getNumberOfSegments()):
                    self.segments[i].addProbability(1.0 * del_prob )

                #tot_prob = 0.0                
                #for i in range(self.getNumberOfSegments()):
                #    tot_prob += self.segments[i].getProbability()
                #print(str(tot_prob))                   
                   
                return
                
            else:
                print('Error: Invalid merge mode')

        # Not enough bins -> split
        if len(self.segments) < self.target_number_of_segments:
            for c in range(self.target_number_of_segments - len(self.segments)):
                split_index = 0
                split_probabilities = []
                for segment in self.segments:
                    split_probabilities.append(segment.getProbability()/prob_tot)
                rand = rnd.random()
                cumulated_probability = 0.0
                for index in range(len(split_probabilities)):
                    cumulated_probability += split_probabilities[index]
                    if cumulated_probability >= rand:
                        split_index = index
                        break
                # We have a split index so do splitting now
                split_segment = self.segments[split_index]
                split_prob = split_segment.getProbability()/2.0
                self.segments[split_index].subProbability(split_prob)
                __segment = Segment(probability         = split_prob,
                                    parent_iteration_id = split_segment.getParentIterationId(),
                                    parent_bin_id       = split_segment.getParentBinId(),
                                    parent_segment_id   = split_segment.getParentSegmentId(),
                                    iteration_id        = split_segment.getIterationId(),
                                    bin_id              = split_segment.getBinId(),
                                    segment_id          = len(self.segments))
                self.__addSegment(__segment)
            return






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
        return segment.getId()

    def __fixSegmentIds(self):
        """
        After deletion of trajectories some segment_ids need to be fixed.
        Call this whenever a segment is deleted from self.segments. This
        should only happen during resampling
        """
        for index in range(len(self.segments)):
            self.segments[index].setSegmentId(index)
    
    def isStartStateBin(self):
        """
        @return Is this bin a start state bin
        """
        return self.is_start_bin
    
    def isEndStateBin(self):
        """
        @return Is this bin a end state bin
        """
        return self.is_end_bin
    
    def getCoordinateIds(self):
        """
        @return list of coordinate ids
        """
        return self.coordinate_ids
    
    def getReferenceNameString(self):
        """
        @return bin reference as a string
        """
        # clumsy construct to use segments own naming method and format
        segment = Segment(probability=0,
                          parent_iteration_id = 0, 
                          parent_bin_id=0, 
                          parent_segment_id=0,
                          iteration_id = self.reference_iteration_id,
                          bin_id       = self.reference_bin_id,
                          segment_id   = self.reference_segment_id)
        return segment.getNameString()
        
    def getNameString(self):
        """
        @return bin index as a string
        """
        name_string = str(self.bin_id).zfill(5)
        return name_string

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
        
    def getInitialProbability(self):
        """
        Returns the initial (before resampling and recycling )
        cumulative probability of all bin trajectories
        """
        probability = 0.0
        for segment in self.initial_segments:
            probability += segment.getProbability()
        return probability
    
    def getNumberOfSegments(self):
        """
        @return Current number of segments
        """
        return len(self.segments)

    def getNumberOfInitialSegments(self):
        """
        @return Number of initial segments
        """
        return len(self.initial_segments)
        
    def getNumberOfPropagatedSegments(self):
        """
        @return Current number of segments if converged is false.
        """
        #if self.isConverged() == False:
        return len(self.segments)
        #else:
        #    return 0
    
    def getTargetNumberOfSegments(self):
        """
        @return Target number of segments
        """
        return self.target_number_of_segments
    
    def __testBinState(self, STATES):
        """
        tests if coordinate ids of bin are in STATES (e.g. start or end state)
        @return boolean
        """
        if type(STATES) == bool:
            return STATES
        
        for state in STATES:
            # compare element wise
            if self.getCoordinateIds() == state:
                return True
        return False
    
    def __iter__(self):
        """
        Defines the class as iterable.
        """
        self._iter_index = -1
        return self

    def next(self):
        """
        Returns the next element of the array self.segments for python 2.6
        """
        self._iter_index += 1
        if self._iter_index >= len(self.segments):
            raise StopIteration
        else:
            return self.segments[self._iter_index]
            
    def __next__(self):
        """
        Returns the next element of the array self.segments
        """
        self._iter_index += 1
        if self._iter_index >= len(self.segments):
            raise StopIteration
        else:
            return self.segments[self._iter_index]

    def backupInitialSegments(self):
        """
        Hard copy the segments list to the initial_segments list
        This function should be called directly after assigning all
        segments to all bins after MD.
        """
        self.initial_segments = copy.deepcopy(self.segments)
        