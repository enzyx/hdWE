#!/usr/bin/python2
import argparse
import sys
import numpy as np
from lib.logger import Logger
import lib.functions_ana_general as functions_ana_general

#####################
###### classes ######
#####################

#####################
##### functions #####
#####################

def assignBinsToStates(iteration, state_A, state_B):
    """
    @return a list of state assignments, sorted according to bin ids. e.g. [A, B, A, '0', A]
    """
    state_list = ['0'] * iteration.getNumberOfBins()
    boundaries = iteration.getBoundaries()
    for this_bin in iteration:
        state_per_dimension = ['0'] * len(boundaries)
        for dimension, coordinate_id in enumerate(this_bin.getCoordinateIds()):
            this_bin_boundaries = []
            if coordinate_id == 0:
                this_bin_boundaries = [0.0, boundaries[dimension][0]]
            elif coordinate_id == len(boundaries[dimension]):
                this_bin_boundaries = [boundaries[dimension][-1], 'inf']
            else:
                this_bin_boundaries = boundaries[dimension][coordinate_id-1:coordinate_id+1]

            if this_bin_boundaries[0] < state_A[dimension][1] and \
               this_bin_boundaries[1] >=  state_A[dimension][0]:
                state_per_dimension[dimension] = 'A'
            elif this_bin_boundaries[0] < state_B[dimension][1] and \
                 this_bin_boundaries[1] >=  state_B[dimension][0]:
                state_per_dimension[dimension] = 'B'                
            
        # collect all state assignments
        this_state = state_per_dimension[0]
        for state in state_per_dimension[1:]:
            if state != this_state:
                this_state = '0'
                break        
        state_list[this_bin.getId()] = this_state
            
    return state_list

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=
    'Bare Model to load iterations. ')
parser.add_argument('-l', '--log', type=str, dest="logdir",
                    required=True, default="hdWE-log", metavar="DIR",
                    help="The logdir to load.")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0, metavar='INT',
                    help="First iteration to use.")
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1, metavar='INT',
                    help="Last iteration to to use.")
parser.add_argument('--state-A', dest="state_A",
                    required=False, type=float, nargs=2, 
                    help="Boundaries of the start state for rate calculation.")  
parser.add_argument('--state-B', dest="state_B",
                    required=False, type=float, nargs=2, 
                    help="Boundaries of the end state for rate calculation.")
parser.add_argument('--output', dest="output_file", 
                    required=False, type=str, default='ana_trace_flux.dat',
                    help="output filename for A and B flux properties")  


# Initialize
args = parser.parse_args()
first_iteration = args.first_iteration
last_iteration  = args.last_iteration
state_A         = [np.sort(args.state_A)] 
state_B         = [np.sort(args.state_B)]

################################

# Get the iterations and parameters
logger = Logger(args.logdir)
if last_iteration < 0:
    last_iteration = logger.getLastIterationId()

current_iteration = logger.loadIteration(first_iteration)
state_list = assignBinsToStates(current_iteration, state_A, state_B)
N = current_iteration.getNumberOfSegments()

# assign initial probabilities
for this_bin in current_iteration:
    for this_segment in this_bin:
        if state_list[this_bin.getId()] == 'A':
            this_segment.setProbability(np.array([1.0/N,0.0]))
        elif state_list[this_bin.getId()] == 'B':
            this_segment.setProbability(np.array([0.0,1.0/N]))
        else:
            this_segment.setProbability(np.array([0.5/N,0.5/N]))
            

flux_into_A         = []
flux_into_B         = []
probability_state_A = []
probability_state_B = []
probability_from_A  = []
probability_from_B  = []

# Iteration Loop
for i in range(first_iteration + 1, last_iteration + 1):
    sys.stdout.write('Iteration {:05d}\r'.format(i))
    sys.stdout.flush()
    previous_iteration = current_iteration
    current_iteration = logger.loadIteration(i)
    state_list = assignBinsToStates(current_iteration, state_A, state_B)
    
    # Initialize data
    flux_into_A_iter         = 0.0 
    flux_into_B_iter         = 0.0
    probability_state_A_iter = 0.0 
    probability_state_B_iter = 0.0 
    probability_from_A_iter  = 0.0 
    probability_from_B_iter  = 0.0  
       
    for this_bin in current_iteration:
        # assign new probability to initial segments
        for this_initial_segment in this_bin.initial_segments:
            new_probability = previous_iteration\
                                .bins[this_initial_segment.getParentBinId()]\
                                .segments[this_initial_segment.getParentSegmentId()].getProbability()
            this_initial_segment.setProbability(new_probability)

        # No split and merge
        if this_bin.getNumberOfSegments() == this_bin.getNumberOfInitialSegments():
            for this_segment in this_bin:
                this_segment.setProbability(this_bin.initial_segments[this_segment.getId()].getProbability())
        
        # Split
        elif this_bin.getNumberOfSegments() > this_bin.getNumberOfInitialSegments():
            split_dict = {}
            for this_segment in this_bin:
                try:
                    split_dict[this_segment.getParentNameString()] += 1
                except KeyError:
                    split_dict[this_segment.getParentNameString()] = 1
            
            for this_initial_segment in this_bin.initial_segments:
                split_dict[this_initial_segment.getParentNameString()] = this_initial_segment.getProbability() / \
                                                split_dict[this_initial_segment.getParentNameString()]
        
            for this_segment in this_bin:
                this_segment.setProbability(split_dict[this_segment.getParentNameString()])
        
        # Merge
        
        elif this_bin.getNumberOfSegments() < this_bin.getNumberOfInitialSegments():
            merge_dict = {}
            kept_prob = 0.0
            merged_probability = 0.0
            # Collect all initial probabilities in a dictionary
            for this_initial_segment in this_bin.initial_segments:
                merge_dict[this_initial_segment.getParentNameString()] = this_initial_segment.getProbability()
            # Assign initial probs to kept segments and set to 0.0 if consumed
            for this_segment in this_bin:
                this_segment.setProbability(merge_dict[this_segment.getParentNameString()])
                kept_prob += this_segment.getProbability()
                merge_dict[this_segment.getParentNameString()] = 0.0
            # Collect the remaining probability
            for parent_string_id in merge_dict:
                merged_probability += merge_dict[parent_string_id]

            merged_probability =  merged_probability/ this_bin.getNumberOfSegments()
            # assign merge probability

            for this_segment in this_bin:
                this_segment.setProbability( this_segment.getProbability() + merged_probability )
                
        # STATES
        if this_bin.getNumberOfSegments() > 0:
            if state_list[this_bin.getId()] == 'A':
                probability_state_A_iter += sum(this_bin.getProbability())
                for this_segment in this_bin:
                    probability = this_segment.getProbability()
                    flux_into_A_iter  += probability[1]
                    this_segment.setProbability(np.array([sum(probability), 0.0])) 
            elif state_list[this_bin.getId()] == 'B':
                probability_state_B_iter += sum(this_bin.getProbability())
                for this_segment in this_bin:
                    probability = this_segment.getProbability()
                    flux_into_B_iter  += probability[0]
                    this_segment.setProbability(np.array([0.0, sum(probability)]))  
                     
    probability_from_A_iter = current_iteration.getProbability()[0]
    probability_from_B_iter = current_iteration.getProbability()[1]
                    
    flux_into_A.append(flux_into_A_iter)         
    flux_into_B.append(flux_into_B_iter) 
    probability_state_A.append(probability_state_A_iter)         
    probability_state_B.append(probability_state_B_iter)   
    probability_from_A.append(probability_from_A_iter)         
    probability_from_B.append(probability_from_B_iter)  


##########################
######### OUTPUT #########
##########################
fout = open(args.output_file, 'w')
fout.write("# F->A           P_A            P->A            F->B            P_B             P->B \n")
# loop over one of the properties and write to file
for i in range(len(flux_into_A)):
    fout.write("{:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}   {:8.7e}\n".format(
                         flux_into_A[i], 
                         probability_state_A[i], 
                         probability_from_A[i], 
                         flux_into_B[i], 
                         probability_state_B[i], 
                         probability_from_B[i]))

fout.close()



print "A --> B"
print   np.mean(flux_into_B[200:])    
print   np.mean(probability_state_A[200:])
print   np.mean(probability_from_A[200:])  
print 'rate:'
print   np.mean(flux_into_B[200:])  / np.mean(probability_state_A[200:])
print '1/mfpt:'
print   np.mean(flux_into_B[200:])  / np.mean(probability_from_A[200:])   


print "B --> A"
print   np.mean(flux_into_A[200:])    
print   np.mean(probability_state_B[200:])
print   np.mean(probability_from_B[200:])  
print 'rate:'
print   np.mean(flux_into_A[200:])  / np.mean(probability_state_B[200:])
print '1/mfpt:'
print   np.mean(flux_into_A[200:])  / np.mean(probability_from_B[200:])    
