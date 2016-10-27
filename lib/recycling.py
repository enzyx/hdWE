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
"""
Steady state code for moving probability from end state bin
to start state bin
"""

def recycleProbability(iteration):
    """
    
    """
    calculateProbabilityFlow(iteration)
    
    if iteration.getProbabilityFlow() > 0.0:
        # Check if segments exist in start state bins
        start_segments = 0
        for start_bin in iteration.getStartStateBins():
            start_segments += start_bin.getNumberOfSegments()
        
        # No segments exist in start state bins, spawn a segment
        # It is assumed that the starting structure lies in the starting state, 
        # then the corresponding bin id is 0
        if start_segments == 0:
            iteration.bins[0].generateSegment(probability         = iteration.getProbabilityFlow(),
                                              parent_iteration_id = 0,
                                              parent_bin_id       = 0,
                                              parent_segment_id   = 0)           
        
        # Segments exist: Scale all segment probabilities accordingly
        else:
            # get normalization factor
            start_state_probability = 0.0
            for start_bin in iteration.getStartStateBins():
                start_state_probability += start_bin.getProbability()
    
            scale_factor = 1.0 + iteration.getProbabilityFlow() / start_state_probability
            
            # scale probabilities of START_STATE segments
            for start_bin in iteration.getStartStateBins():
                for this_segment in start_bin:
                        this_segment.multProbability(scale_factor)   
                        
    # Empty all end state bins
    emptyEndStateBins(iteration)

def calculateProbabilityFlow(iteration):
    """
    Sum up the probabilities of all segments in end state bins
    and save the value in the probability_flow variable of the iteration
    """
    probability_flow = 0.0
    for end_state_bin in iteration.getEndStateBins():
        for end_state_segment in end_state_bin:
            probability_flow += end_state_segment.getProbability()
    
    # reassign flown probability
    iteration.setProbabilityFlow(probability_flow)
    
def emptyEndStateBins(iteration):
    """
    Delete all segments in the end state bins (their probability 
    is now assigned to the start state bins)
    """
    for end_state_bin in iteration.getEndStateBins():
        end_state_bin.segments = []                  
