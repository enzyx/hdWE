from lib.segment import Segment
import numpy
import random as rnd

class Resampling(object):
    """
    Contain all resampling modes and 
    """
    def __init__(self, 
                 md_module, 
                 resampling_mode, 
                 closest_merge_threshold = 0.0,
                 split_forward_coordinate_id = 0,
                 split_forward_number_of_children = 1):
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
        
        # Split forward mode
        self.SPLIT_FORWARD_COORDINATE_ID      = split_forward_coordinate_id
        self.SPLIT_FORWARD_NUMBER_OF_CHILDREN = split_forward_number_of_children

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

        elif self.RESAMPLING_MODE == 'weighted':
            self.resampleWeighted()
        
        elif self.RESAMPLING_MODE == 'random':
            self.resampleRandom()
        
        elif self.RESAMPLING_MODE == 'no-merge':
            self.resampleNoMerge()
    
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
            if this_bin.getNumberOfSegments() < this_bin.target_number_of_segments:
                if this_bin.getId() == 0:
                    self.splitWeighted(this_bin)
                else:
                    self.splitForward(this_bin)
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
            if this_bin.getNumberOfSegments() > this_bin.target_number_of_segments:
                self.mergeWeighted(this_bin)
                continue
            
            # Not enough segments -> split
            if this_bin.getNumberOfSegments() < this_bin.target_number_of_segments:
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
            if this_bin.getNumberOfSegments() > this_bin.target_number_of_segments:
                self.mergeRandom(this_bin)
                continue
            
            # Not enough segments -> split
            if this_bin.getNumberOfSegments() < this_bin.target_number_of_segments:
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
            if this_bin.getNumberOfSegments() < this_bin.target_number_of_segments:
                self.splitWeighted(this_bin)
                continue

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
            if self.iteration.bins[parent_bin_id].getCoordinateIds()[self.SPLIT_FORWARD_COORDINATE_ID] > \
                                         this_bin.getCoordinateIds()[self.SPLIT_FORWARD_COORDINATE_ID]:
                # split this segment
                split_indices[this_segment.getId()] = self.SPLIT_FORWARD_NUMBER_OF_CHILDREN
        
        # Reduce the list
        reduced_split_indices = []
        for segment_id, n_children in enumerate(split_indices):
            if n_children > 0:
                reduced_split_indices.append([segment_id, n_children])
        
        while len(reduced_split_indices) > 0 and this_bin.getNumberOfSegments < this_bin.target_number_of_segments:
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
            for this_segment in this_bin.segments:
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