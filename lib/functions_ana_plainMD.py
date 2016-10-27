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
import numpy
import lib.constants as constants
from math import log

def PMF(data, N_bins):
    """
    Returns the PMF and probability distribution corresponding to data.
    The energy unit is according to constants.py
    """
    # Calculate the weighted histogram and PMF     
    # Setup variables
    hist_min =  min(data)
    hist_max =  max(data)
    
    dcoord   =  1.0 * (hist_max - hist_min ) / N_bins
    hist     =  numpy.zeros([N_bins,3], float)
    # Sort coords into histogramm
    for i in range(0,len(data)):
        index       = int( (data[i] - hist_min) / dcoord )
        #maximum coord entry shall be included in last bin:
        if index==N_bins:
            index = index - 1
        hist[index,2] = hist[index,2] + 1
    # Assign the bin positions and calculate free energy:
    for i in range(0, N_bins):
        hist[i,0] = hist_min + (i + 0.5) * dcoord
        if hist[i,2] > 0:
            hist[i,1]  = - constants.kT * log(hist[i,2])
        else:
            hist[i,1]  = 'Inf'
    
    # Shift minimum to zero        
    pmf_min = min(hist[:,1])
    for i in range(0, N_bins):
        hist[i,1] = hist[i,1] - pmf_min
        
    return hist

def inState(coordinate, state):
    """
    Returns True if coordinate lies within the state boundaries [lower_boundary, upper_boundary),
    otherwise False.
    state has to be given as:
    state[0] = lower_boundary
    state[1] = upper_boundary
    """
    if (coordinate > state[0]) and (coordinate < state[1]):
        return True
    else:
        return False
        
def rate(N_transitions_into_target_state, total_residence_time_in_start_state):
    """
    Returns the rate: number of transitions from start to target state
    devided by the total residence time in the start state
    """
    return float(N_transitions_into_target_state) / float(total_residence_time_in_start_state)


def transitions_from_coordinates(coordinates, start_state, end_state):
    """
    @ return list with length transitions 
    that holds the residence times before this transition
    """
    transitions   = [] # list of with residence times prior to each transition
    in_transition = False
    
    for coord in coordinates:
        if in_transition == False:
            # if in start state, start with transition time counting
            if inState(coord, start_state):
                in_transition  = True
                residence_time = 1
                
        if in_transition == True:
            # if in transition, save transition time if arrived in end state,
            # otherwise add frame to transition time 
            if inState(coord, end_state):
                transitions.append(residence_time)
                in_transition    = False
            if inState(coord, start_state):
                residence_time += 1

    if len(transitions) == 0:
        print(' No transitions found.')
      
    return transitions

def mfpt_from_coordinates(coordinates, start_state, end_state):
    first_passage_times    = []
    in_transition          = False
    first_passage_time_tmp = 0
    
    for coord in coordinates:
        if in_transition == False:
            # if in start state, start with transition time counting
            if inState(coord, start_state):
                in_transition          = True
                first_passage_time_tmp = 1
                
        if in_transition == True:
            # if in transition, save transition time if arrived in end state,
            # otherwise add frame to transition time 
            # and keep track of new possible transition starting points
            if inState(coord, end_state):
                first_passage_times.append(first_passage_time_tmp)
                in_transition    = False
            else:
                first_passage_time_tmp += 1

    if len(first_passage_times) == 0:
        print(' No transitions found.')
      
    return first_passage_times

  
def residence_times_from_coordinates(coordinates, state):

    residence_times        = []
    residence_time_tmp     = 0
    already_in_state       = False
    
    for coord in coordinates:
        if already_in_state == False:
            if inState(coord, state):
                residence_time_tmp = 1
                already_in_state   = True
        else:
            if inState(coord, state):
                residence_time_tmp += 1
            else:
                already_in_state   = False
                residence_times.append(residence_time_tmp)
           
    return residence_times
    
def distribution_from_times(transition_times, N_bins):
    distr_min = min(transition_times)
    distr_max = max(transition_times)
    d         =  1.0 * (distr_max - distr_min ) / N_bins
    distr     =  numpy.zeros([N_bins,2], float)
    #Sort coords into histogramm
    for i in range(0,len(transition_times)):
        index       = int( (transition_times[i] - distr_min) / d )
        #maximum coord entry shall not be in an extra bin:
        if index==N_bins:
            index = index - 1
        distr[index,1] += 1
    #Assign the bin positions and normalize:
    for i in range(0,N_bins):
        distr[i,0] = distr_min + i * d
        distr[i,1] /= len(distr)
   
    return distr
    
def F_of_state_from_hist(hist, state):
    p = 0.0
    for i in range(len(hist)):
        if (hist[i,0] >= state[0]) and (hist[i,0] < state[1]):
            p += hist[i,2]
    return - constants.kT * log(p)
