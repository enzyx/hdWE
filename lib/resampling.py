#
# This file is part of hdWE. 
# Copyright (C) 2016 Manuel Luitz <manuel.luitz@tum.de>
# Copyright (C) 2016 Rainer Bomblies <r.bomblies@tum.de>
# Copyright (C) 2016 Fabian Zeller
#
# hdWE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# hdWE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with hdWE. If not, see <http://www.gnu.org/licenses/>.
# 
from lib.segment import Segment
import numpy
import random as rnd
from resampling_events import Merge,Split

class Resampling(object):
    """
    Contain all resampling modes and 
    """
    def __init__(self, 
                 md_module, 
                 resampling_mode, 
                 closest_merge_threshold = 0.0,
                 primary_coordinate  = 0,
                 split_forward_number_of_children = 1,
                 split_region = [0, 9e99],
                 front_interval = 9e99,
                 west_weight_split_threshold=2.0,
                 west_weight_merge_cutoff=1.0,
                 west_do_adjust_counts=True):
        """
        Keep track of some resampling variables
        """
        # a local copy to the iteration object (updated when 
        # self.resample() is called)
        self.iteration                        = None
        self.md_module                        = md_module
        
        # Merge mode
        self.RESAMPLING_MODE                  = resampling_mode
        
        # Closest mode
        self.CLOSEST_MERGE_THRESHOLD          = closest_merge_threshold
        
        # Split forward (front) mode
        self.PRIMARY_COORDINATE               = primary_coordinate
        self.SPLIT_FORWARD_NUMBER_OF_CHILDREN = split_forward_number_of_children
        self.FRONT_INTERVAL                   = front_interval
        
        # Split region
        self.SPLIT_REGION                     = split_region
        
        # westpa constants, default values are taken 
        # from the WESTPA code. 
        self.WESTPA_WEIGHT_SPLIT_THRESHOLD    = west_weight_split_threshold
        self.WESTPA_DO_ADJUST_COUNTS          = west_do_adjust_counts
        self.WESTPA_WEIGHT_MERGE_CUTOFF       = west_weight_merge_cutoff

    def resample(self, iteration):
        """
        Split or Merge segments to generate the target number of segments
        """
        # Update the link to iteration
        self.iteration = iteration
        
        if self.RESAMPLING_MODE == 'closest':
            bins_rmsds = self.md_module.calcBinsRmsds(self.iteration)
            self.resampleClosest(bins_rmsds)
                            
        elif self.RESAMPLING_MODE == 'split-forward':
            self.resampleSplitForward()
            
        elif self.RESAMPLING_MODE == 'split-forward-front':
            front_bin_coord_id = 9e99
            for this_bin in iteration.bins:
                if this_bin.getNumberOfSegments() > 0:
                    if this_bin.getCoordinateIds()[self.PRIMARY_COORDINATE] < front_bin_coord_id:
                        front_bin_coord_id = this_bin.getCoordinateIds()[self.PRIMARY_COORDINATE]
            self.resampleSplitForwardFront(front_bin_coord_id)
        
        elif self.RESAMPLING_MODE == 'split-region':
            self.resampleSplitRegion()
    
        elif self.RESAMPLING_MODE == 'weighted':
            self.resampleWeighted()
        
        elif self.RESAMPLING_MODE == 'random':
            self.resampleRandom()
        
        elif self.RESAMPLING_MODE == 'no-merge':
            self.resampleNoMerge()
        
        elif self.RESAMPLING_MODE == 'westpa':
            self.resampleWestpa()
    
    ##########################################
    #    Implementation of resampling modes  #
    ##########################################
    def resampleClosest(self, bins_rmsds):
        """
        This mode merges segments according to their conformations RMSD.
        The closest conformations are merged together until the target number of 
        segments is reached.  
        """
        for this_bin in self.iteration:
            rmsds = bins_rmsds[this_bin.getId()]
            
            # Is in sample region or is bin empty?
            if len(this_bin.segments) == 0 or this_bin.sample_region == False:
                continue
            
            # To many segments -> merge
            if this_bin.getNumberOfSegments() > this_bin.target_number_of_segments:
                self.mergeClosest(this_bin, rmsds)
                continue
            
            # Not enough segments -> split
            if this_bin.getNumberOfSegments() < this_bin.target_number_of_segments:
                self.splitWeighted(this_bin)
                continue

    def resampleSplitForward(self):
        """
        This resampling mode only splits segments when they move forward along 
        a given coordinate.
        """
        for this_bin in self.iteration:
            # Is in sample region or is bin empty?
            if len(this_bin.segments) == 0 or this_bin.sample_region == False:
                continue
                        
            # Not enough segments -> split
            if this_bin.getNumberOfSegments() < this_bin.getTargetNumberOfSegments():
                if this_bin.getId() == 0:
                    self.splitWeighted(this_bin)
                else:
                    self.splitForward(this_bin)
                continue
            
    def resampleSplitForwardFront(self, front_bin_coord_id):
        """
        This resampling mode only splits segments when they move forward along 
        a given coordinate. 
        Segments in bins which are not in the front are merged until the corresponding bin ends up with
        a max of 1 segment.
        Front bins are defined by a range from the bin with the lowest coordinate id.        
        """
        for this_bin in self.iteration:
            # Is in sample region or is bin empty?
            if len(this_bin.segments) == 0 and this_bin.getId() > 0:
                continue
            if this_bin.sample_region == False:
                continue
             
            # Bin is in front region           
            if this_bin.getCoordinateIds()[self.PRIMARY_COORDINATE] < front_bin_coord_id + self.FRONT_INTERVAL:
                # Not enough segments -> split
                if this_bin.getNumberOfSegments() < this_bin.getTargetNumberOfSegments():
                    if this_bin.getId() == 0:
                        self.splitFromEnsemble(this_bin)
                    else:
                        self.splitForward(this_bin)
                    continue
            # Bin is not in front region
            else:
                if this_bin.getId() == 0 and this_bin.getNumberOfSegments() == 0:
                    print('Bin is not in front region!')
                    self.splitFromEnsemble(this_bin, True)
                else:
                    self.mergeKeepRandomSegment(this_bin)
    
    def resampleSplitRegion(self):
        """
        This resampling mode only splits segments in a given region of bins.
        """
        for this_bin in self.iteration:
            # Is in sample region or is bin empty?
            if len(this_bin.segments) == 0 or this_bin.sample_region == False:
                continue
                        
            # Not enough segments -> split
            if this_bin.getNumberOfSegments() < this_bin.getTargetNumberOfSegments():
                if self.SPLIT_REGION[0] <= this_bin.getCoordinateIds()[self.PRIMARY_COORDINATE] <= self.SPLIT_REGION[1]:
                    self.splitWeighted(this_bin)
                elif this_bin.getId() == 0:
                    self.splitWeighted(this_bin)
                continue
    
    def resampleWeighted(self):
        """
        This mode merges segments by its weighted probabilities and splits 
        by the same method.
        """
        for this_bin in self.iteration:
            # Is in sample region or is bin empty?
            if len(this_bin.segments) == 0 or this_bin.sample_region == False:
                continue
            
            # To many segments -> merge
            if this_bin.getNumberOfSegments() > this_bin.getTargetNumberOfSegments():
                self.mergeWeighted(this_bin)
                continue
            
            # Not enough segments -> split
            if this_bin.getNumberOfSegments() < this_bin.getTargetNumberOfSegments():
                self.splitWeighted(this_bin)
                continue

    def resampleRandom(self):
        """
        This mode merges segments randomly and splits segments weighted. 
        """
        for this_bin in self.iteration:
            # Is in sample region or is bin empty?
            if len(this_bin.segments) == 0 or this_bin.sample_region == False:
                continue
            
            # To many segments -> merge
            if this_bin.getNumberOfSegments() > this_bin.getTargetNumberOfSegments():
                self.mergeRandom(this_bin)
                continue
            
            # Not enough segments -> split
            if this_bin.getNumberOfSegments() < this_bin.getTargetNumberOfSegments():
                self.splitWeighted(this_bin)
                continue

    def resampleNoMerge(self):
        """
        Only split (weighted) segments but do not merge
        """
        for this_bin in self.iteration:
            # Is in sample region or is bin empty?
            if len(this_bin.segments) == 0 or this_bin.sample_region == False:
                continue
            
            # No merging
            
            # Not enough segments -> split
            if this_bin.getNumberOfSegments() < this_bin.getTargetNumberOfSegments():
                self.splitWeighted(this_bin)
                continue
    
    def resampleWestpa(self):
        '''
        Resampling mode as implemented in the WESTPA code.
        Comments are sometimes take from the original implementation.
        '''       
        # Regardless of current particle count, always split overweight 
        # segments and merge underweight particles.
        # Then and only then adjust for correct particle count
        for this_bin in self.iteration:
            if len(this_bin.segments) == 0 or this_bin.sample_region == False:
                continue
            
            self._westpaSplitByWeight(this_bin)
            self._westpaMergeByWeight(this_bin)
            if self.WESTPA_DO_ADJUST_COUNTS:
                self._westpaAdjustCount(this_bin)

    def _westpaSplitByWeight(self, this_bin):
        '''Split overweight particles'''
        weights      = this_bin.getProbability()
        target_count = this_bin.getTargetNumberOfSegments()
        ideal_weight = weights / target_count
        
        to_split     = []
        for segment in this_bin:
            if segment.getProbability() > self.WESTPA_WEIGHT_SPLIT_THRESHOLD * ideal_weight:
                to_split.append(segment)
        
        for segment in to_split:
            m = int(numpy.ceil(segment.getProbability() / ideal_weight))
            self._westpaSplitWalker(segment, m, this_bin)

    def _westpaSplitWalker(self, segment, m, this_bin):
        '''Split the walker ``segment`` (in ``bin``) into ``m`` walkers
           Segment has therefore m-1 children
        '''
        split_prob = segment.getProbability()/float(m)
        this_bin.segments[segment.getId()].setProbability(split_prob)
        
        for dummy in range(1, m):
            __segment = Segment(probability         = split_prob,
                                parent_iteration_id = segment.getParentIterationId(),
                                parent_bin_id       = segment.getParentBinId(),
                                parent_segment_id   = segment.getParentSegmentId(),
                                iteration_id        = segment.getIterationId(),
                                bin_id              = segment.getBinId(),
                                segment_id          = this_bin.getNumberOfSegments())
            this_bin.addSegment(__segment)
        
        # Save history
        this_bin.resampling_history.append(Split(segment.getId(), m))
    
    def _westpaMergeByWeight(self, this_bin):
        '''Merge underweight particles'''
        weight       = this_bin.getProbability()
        target_count = this_bin.getTargetNumberOfSegments()
        ideal_weight = weight / target_count

        while True:
            # Sorted segments weights
            sorted_segs = sorted(this_bin.segments[:], key=lambda seg: seg.getProbability())
                        
            to_merge = []
            cumul_weight = 0.0
            for segment in sorted_segs:
                cumul_weight += segment.getProbability()
                if cumul_weight <= ideal_weight * self.WESTPA_WEIGHT_MERGE_CUTOFF:
                    to_merge.append(segment)
                else: break
            
            if len(to_merge) < 2:
                return
            
            self._westpaMergeWalkers(to_merge, None, this_bin)
    
    def _westpaMergeWalkers(self, segments, cumul_weight, this_bin):
        """
        Merge the given ``segments`` in ``bin``, previously sorted by weight, 
        into one conglomerate segment. ``cumul_weight`` is the cumulative sum of the weights 
        of the ``segments``; this may be None to calculate here.
        """
        if cumul_weight is None:
            cumul_weight = numpy.add.accumulate([segment.getProbability() for segment in segments])
        # Total probability
        tot_weight = cumul_weight[-1]
        # Get index of the surviving segment
        iparent     = numpy.digitize((rnd.uniform(0, tot_weight),), cumul_weight)[0]
        psegment    = segments[iparent]
        del segments[iparent]
        iextinction = [segment.getId() for segment in segments]
        
        # Save history original segment ids
        this_bin.resampling_history.append(Merge(psegment.getId(), iextinction))
        
        # Assign all probability to surviving segment
        psegment.setProbability(tot_weight)
        
        iextinction.sort(reverse=True)
        for index in iextinction:
            del this_bin.segments[index]
        this_bin.fixSegmentIds()

        
    
    def _westpaAdjustCount(self, this_bin):
        # Split
        while this_bin.getNumberOfSegments() < this_bin.getTargetNumberOfSegments():
            # Always split the highest probability walker into two
            split_segment = max(this_bin.segments, key=lambda segment: segment.getProbability())
            self._westpaSplitWalker(split_segment, 2, this_bin)
        
        # Merge
        while this_bin.getNumberOfSegments() > this_bin.getTargetNumberOfSegments():
            # Create a copy of references
            sorted_segs = sorted(this_bin.segments[:], key=lambda seg: seg.getProbability())
            self._westpaMergeWalkers(sorted_segs[:2], None, this_bin)
            
        

######################################################
#   PARALLELIZATION CODE FOR RESAMPLING. 
#   CAN BE USED AGAIN IF DESIRED.    
#       # Parallel
#     if NUMBER_OF_THREADS > 1:
#         thread_container = ThreadContainer()
#         for this_bin in iterations[-1]:
#             thread_container.appendJob(threading.Thread(target = this_bin.resampleSegments,  
#                                                         args= (MERGE_MODE, MERGE_THRESHOLD)))
#             if thread_container.getNumberOfJobs() >= NUMBER_OF_THREADS:
#                 thread_container.runJobs()
#         # Run remaining jobs
#         thread_container.runJobs()
#     # Serial
#     else:
#         for this_bin in iterations[-1]:
#             this_bin.resampleSegments(MERGE_MODE, 
#                                       MERGE_THRESHOLD)
#  
######################################################
   
    ########################
    #       Split modes    #
    ########################
    def splitWeighted(self, this_bin):
        """
        Select a segment randomly according to its probability and
        split it. Repeat the process until we reach the 
        target number of segments. Segments with high probability
        have a higher chance to be split than segments with lower 
        probability.
        """
        prob_tot = this_bin.getProbability()
        split_indices = [0] * this_bin.getNumberOfSegments()
        split_probabilities = []
        for segment in this_bin.segments:
            split_probabilities.append(segment.getProbability()/prob_tot)

        for dummy in range(this_bin.target_number_of_segments - len(this_bin.segments)):
            rand = rnd.random()
            cumulated_probability = 0.0
            for index in range(len(split_probabilities)):
                cumulated_probability += split_probabilities[index]
                if cumulated_probability >= rand:
                    split_indices[index] += 1
                    break
        
        # We have a list how often each segment is split
        for segment_id, number_of_children in enumerate(split_indices):
            if number_of_children == 0:
                continue
            split_segment = this_bin.segments[segment_id]
            split_prob = split_segment.getProbability()/float(number_of_children + 1)
            
            this_bin.segments[segment_id].setProbability(split_prob)
            for dummy in range(number_of_children):
                __segment = Segment(probability         = split_prob,
                                    parent_iteration_id = split_segment.getParentIterationId(),
                                    parent_bin_id       = split_segment.getParentBinId(),
                                    parent_segment_id   = split_segment.getParentSegmentId(),
                                    iteration_id        = split_segment.getIterationId(),
                                    bin_id              = split_segment.getBinId(),
                                    segment_id          = this_bin.getNumberOfSegments())
                this_bin.addSegment(__segment)

    def splitForward(self, this_bin):
        """
        Split segments which moved along the given coordinate during the
        last iteration.
        """
        split_indices = [0] * this_bin.getNumberOfSegments()
        for this_segment in this_bin.segments:
            parent_bin_id = this_segment.getParentBinId()
            # Check whether this segment came from an earlier bin:
            if self.iteration.bins[parent_bin_id].getCoordinateIds()[self.PRIMARY_COORDINATE] > \
                                         this_bin.getCoordinateIds()[self.PRIMARY_COORDINATE]:
                # split this segment
                split_indices[this_segment.getId()] = self.SPLIT_FORWARD_NUMBER_OF_CHILDREN
        
        # Reduce the list
        reduced_split_indices = []
        for segment_id, n_children in enumerate(split_indices):
            if n_children > 0:
                reduced_split_indices.append([segment_id, n_children])
        
        while len(reduced_split_indices) > 0 and this_bin.getNumberOfSegments() < this_bin.getTargetNumberOfSegments():
            rand_index = rnd.randint(0, len(reduced_split_indices)-1)
            segment_id = reduced_split_indices[rand_index][0]
            n_children = reduced_split_indices[rand_index][1]

            split_segment = this_bin.segments[segment_id]
            split_prob    = split_segment.getProbability()/float(n_children + 1)
            
            this_bin.segments[segment_id].setProbability(split_prob)
            for dummy in range(n_children):
                __segment = Segment(probability         = split_prob,
                                    parent_iteration_id = split_segment.getParentIterationId(),
                                    parent_bin_id       = split_segment.getParentBinId(),
                                    parent_segment_id   = split_segment.getParentSegmentId(),
                                    iteration_id        = split_segment.getIterationId(),
                                    bin_id              = split_segment.getBinId(),
                                    segment_id          = this_bin.getNumberOfSegments())
                this_bin.addSegment(__segment)
            
            del(reduced_split_indices[rand_index])
            
    def splitFromEnsemble(self, this_bin, generate_only_one_segment = False):
        """
        randomly add segments to a bin (usually starting bin with bin_id = 0)
        from the starting structures
        """
        target_number_of_segments = this_bin.getTargetNumberOfSegments()
        if generate_only_one_segment == True:
            target_number_of_segments = 1
        for dummy in range(target_number_of_segments - len(this_bin.segments)):
            random_index = rnd.randint(0, self.iteration.number_starting_structures -1)
            segment = Segment(probability         = 1.0,
                              parent_iteration_id = 0,
                              parent_bin_id       = 0,
                              parent_segment_id   = random_index,
                              iteration_id        = this_bin.getIterationId(),
                              bin_id              = this_bin.getId(),
                              segment_id          = this_bin.getNumberOfSegments())
            this_bin.addSegment(segment)
        

    
    ########################
    #      Merge  modes    #
    ########################
    def mergeClosest(self, this_bin, rmsds):
        """
        Select the pair of segments with the lowest rmsd (closest conformations)
        and eliminate randomly one of both. Repeat until target_number_of_segments
        is reached.
        """
        # deal with the rmsds to the segment itself which is always zero
        for i in range(0, len(rmsds)):
            rmsds[i,i] = 'Inf'
        # check for remaining 0 entries (which are most probably due to RMSD calculation error)
        # Or when to simulations are run with the same random seed (very rare)
        #TODO: remove this message or put it to debugging
        for i in range(0, len(rmsds)):
            for j in range(0, len(rmsds)):
                if rmsds[i,j] == 0.0:
                    print('RMSD is zero, this might be due to an RMSD calculation error')
                    print('or because two MDs were run with the same random seed')
                    print('Here is the RMSD matrix:')
                    print(rmsds)
        for dummy in range(len(this_bin.segments) - this_bin.target_number_of_segments):
            # do not merge segments if the lowest rmsd is above threshold
            if (self.CLOSEST_MERGE_THRESHOLD > 0.0) and (numpy.min(rmsds) > self.CLOSEST_MERGE_THRESHOLD):
                print('     Number of segments in bin {bin:05d} is {n_segs} after merge.'
                      .format(bin = this_bin.getId(), n_segs = this_bin.getNumberOfSegments()) )
                break
            # get indices of the two segments with the lowest RMSD with respect to each other
            # randomly choose one segment of these two that will be merged into the other
            lowest_rmsd_indices  = numpy.unravel_index(numpy.argmin(rmsds), rmsds.shape)
            lowest_rmsd          = rmsds[lowest_rmsd_indices[0], lowest_rmsd_indices[1]]
            target_random_int    = rnd.randint(0,1)
            ext_random_int       = int(not target_random_int)
            target_index = lowest_rmsd_indices[target_random_int]
            ext_index    = lowest_rmsd_indices[ext_random_int]
            # save indices in merge_list
            this_bin.merge_list.append([ext_index, target_index])
            # save the rmsd between the merged structure in merge_rmsd_list
            this_bin.merge_rmsd_list.append(lowest_rmsd)
            # shift the probability of deleted segment to the target segment
            this_bin.segments[target_index].addProbability(this_bin.segments[ext_index].getProbability())
            # delete segments and fix segment ids
            del this_bin.segments[ext_index]
            this_bin.fixSegmentIds()                                        
            # delete corresponding entries from the rmsd matrix
            rmsds = numpy.delete(rmsds, ext_index, 0)
            rmsds = numpy.delete(rmsds, ext_index, 1)
    
    def mergeWeighted(self, this_bin):
        """
        Choose a segment randomly according to probability
        Distribute the probability of the deleted segment
        equally amongst the remaining segments
        """
        extinction_probability = 0.0
        prob_tot = this_bin.getProbability()
        for dummy in range(len(this_bin.segments) - this_bin.target_number_of_segments):
            # Get the extinction index
            ext_index = 0
            inv_weights = []
            inv_weights_tot = 0.0
            for segment in this_bin.segments:
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
                
            # Reassign the extinction probability / N_segments to the remaining segments
            extinction_probability = this_bin.segments[ext_index].getProbability()
            
            this_bin.merge_list.append([ext_index])
            del this_bin.segments[ext_index]
            for this_segment in this_bin.segments:
                this_segment.addProbability(extinction_probability / this_bin.getNumberOfSegments())  
                this_bin.merge_list[-1].append(this_segment.getId())
            # Reorder segment ids after deletion 
            this_bin.fixSegmentIds()
    
    def mergeRandom(self, this_bin):
        """
        Choose a segment randomly. Distribute the 
        probability of the deleted segments
        equally amongst the remaining segments
        """
        extinction_probability = 0.0
        for dummy in range(len(this_bin.segments) - this_bin.target_number_of_segments):
            # Get the extinction index
            ext_index = rnd.randint(0, this_bin.getNumberOfSegments() - 1)
            # Reassign the extinction probability / N_segments to the remaining segments
            extinction_probability = this_bin.segments[ext_index].getProbability()
            this_bin.merge_list.append([ext_index])
            del this_bin.segments[ext_index]
            for this_segment in this_bin.segments:
                this_segment.addProbability(extinction_probability / this_bin.getNumberOfSegments())  
            for this_segment in this_bin.segments:
                this_bin.merge_list[-1].append(this_segment.getId())
        
            # Reorder segment ids after deletion 
            this_bin.fixSegmentIds()
            
    def mergeKeepRandomSegment(self, this_bin):
        """
        Choose a segment randomly. Distribute the 
        probability of the deleted segments
        equally amongst the remaining segments
        """
        extinction_probability = 0.0
        for dummy in range(len(this_bin.segments) - 1):
            # Get the extinction index
            ext_index = rnd.randint(0, this_bin.getNumberOfSegments() - 1)
            # Reassign the extinction probability / N_segments to the remaining segments
            extinction_probability = this_bin.segments[ext_index].getProbability()
            this_bin.merge_list.append([ext_index])
            del this_bin.segments[ext_index]
            for this_segment in this_bin.segments:
                this_segment.addProbability(extinction_probability / this_bin.getNumberOfSegments())  
            for this_segment in this_bin.segments:
                this_bin.merge_list[-1].append(this_segment.getId())
        
            # Reorder segment ids after deletion 
            this_bin.fixSegmentIds()
        
