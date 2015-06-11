#!/usr/bin/env python2
"""
Compare MD and hdWE rama bin properties 
"""
from __future__ import print_function
import numpy
from logger import Logger
import argparse  
import analysis_operations
from numpy.core.test_rational import denominator

# [ {'id':<int>, 'rama_id':'11121', 'p_md':<float>, 'p_hdWE':<float>},
#   {}, 
#   ... ]
# 
bin_data = []

def getIndex(rama_id):
    """
    returns the index of the bin_data element with the given rama_id
    """
    for index, loop_bin_data in enumerate(bin_data):
        if loop_bin_data['rama_id'] == rama_id:
            return index
    return None

def resortMatrix(matrix, permutation_array):
    N = len(matrix)
    sorted_matrix = numpy.zeros((N,N), tuple)
    for i in range(N):
        for j in range(N):
            sorted_matrix[i,j] = matrix[permutation_array[i],permutation_array[j]]
    return sorted_matrix

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    metavar="FILE", required=True,
                    help="The log directory for reading")                  
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0,
                    help="First iteration to use.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Last iteration to use.")
parser.add_argument('-f', '--md-file', type=str, dest="md_rama_bins", 
                    metavar="FILE", default="rama_bins.dat",
                    help="Rama bins file from MD.")
sorter = parser.add_mutually_exclusive_group(required=False)
sorter.add_argument('-m', '--md-sort', dest="md_sort", 
                    default=False, action='store_true',
                    help="Sort bins by MD bin probability.")
sorter.add_argument('-w', '--hdWE-sort', dest="hdWE_sort", 
                    default=False, action='store_true',
                    help="Sort bins by hdWE bin probability.")

# parser.add_argument('-t', '--bin-to-bin-transitions', action="store_true",
#                     dest="bin_to_bin_transitions", 
#                     help="Print a matrix with bin to bin transition events")
# parser.add_argument('-o', '--output', dest="outfile", metavar="OUTPUT",
#                     help="Write the matrix to this file")

args = parser.parse_args()
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration)

# read hdWE data into bin_data
for loop_bin in iterations[-1]:
    bin_data.append({'id':      loop_bin.getId(),
                     'rama_id': loop_bin.getRamaId(),
                     'p_hdWE':  loop_bin.getProbability(),
                     'p_md':    0.0})
    
# read MD data into bin_data

md_rama_frames_file = open(args.md_rama_bins,'r')
nframes = 0
md_rama_frames = []
# Read in the rama bins of the MD in form [[frame_id, rama_id], ... ]
for line in md_rama_frames_file.readlines():
    md_rama_frames.append(line.strip().split())
md_rama_frames_file.close()

# Sort the MD frames into the hdWE bin data structure
# and count the probabilities
for frame in md_rama_frames:
    nframes += 1.0
    frame_id = frame[0]
    frame_rama_id = frame[1]
    frame_bin_index = getIndex(frame_rama_id)
    if frame_bin_index != None:
        bin_data[frame_bin_index]['p_md'] += 1.0
    else:
        bin_data.append({'id':      len(bin_data),
                         'rama_id': frame_rama_id,
                         'p_hdWE':  0.0,
                         'p_md':    1.0})
# Normalize
for loop_bin in bin_data:
    loop_bin['p_md'] /= nframes

# Get hdWE transition and rate matrices 
hdWE_rate_matrix       = analysis_operations.meanRateMatrix(iterations)
hdWE_transition_matrix = analysis_operations.cumulativeTransitionMatrix(iterations)

#########
#M = len(hdWE_transition_matrix)
#hdWE_rate_matrix = numpy.zeros((M,M),float)
#N_segments_bin = iterations[-1].bins[0].getTargetNumberOfSegments()
#tau = 1.0
#N_iter = len(iterations)
#denominator = float(N_segments_bin * tau * N_iter)
#for i in range(M):
#    for j in range(M):
#        hdWE_rate_matrix[i,j] = hdWE_transition_matrix[i,j] / denominator
########

# Get MD transition matrix
N = len(bin_data)
md_transition_matrix = numpy.zeros((N,N), int)
prev_frame = md_rama_frames[0]
for frame in md_rama_frames[1:]:
    i = getIndex(prev_frame[1])
    j = getIndex(frame[1])
    md_transition_matrix[i,j] += 1
    prev_frame = frame

# Get MD rate matrix
md_rate_matrix = numpy.zeros((N,N), float)
for i in range(N):
    T_i = float(sum(md_transition_matrix[i,:]))
    #T_i = float(md_transition_matrix[i,i])
    if T_i > 0.0: 
        for j in range(N):
            md_rate_matrix[i,j] = float(md_transition_matrix[i,j]) / ( T_i )
        
# Combine matrices to a matrix with tuple elements i,j (hdWE_matrix(i,j), md_matrix(i,j))
combined_transition_matrix = numpy.zeros((N,N), tuple)
combined_rate_matrix = numpy.zeros((N,N), tuple)
M = len(hdWE_transition_matrix)
for i in range(N):
    for j in range(N):
        if i >= M or j >= M: 
            combined_transition_matrix[i,j] = (0, md_transition_matrix[i,j])
            combined_rate_matrix[i,j]       = (0, md_rate_matrix[i,j])
        else:
            combined_transition_matrix[i,j] = (hdWE_transition_matrix[i,j], md_transition_matrix[i,j])
            combined_rate_matrix[i,j]       = (hdWE_rate_matrix[i,j], md_rate_matrix[i,j])      
        
# Sort bins
if args.md_sort:
    bin_data = sorted(bin_data, key=lambda e:e['p_md'], reverse=True)
if args.hdWE_sort:
    bin_data = sorted(bin_data, key=lambda e:e['p_hdWE'], reverse=True)

sort_indices = []
for this_bin in bin_data:
    sort_indices.append(this_bin['id'])
    
# Sort combined matrix
combined_transition_matrix = resortMatrix(combined_transition_matrix, sort_indices)
combined_rate_matrix = resortMatrix(combined_rate_matrix, sort_indices)

# Print output
bin_file = 'bin_probabilities.dat'
f = open(bin_file, 'w')
f.write("# rama_id   p_md   p_hdWE   bin_id\n")
for b in sorted(bin_data, key=lambda e:e['p_md'], reverse=True):
    f.write("{0:s} {1:15e} {2:15e} {3:d}\n".format(b["rama_id"], b["p_md"], b["p_hdWE"], b["id"]))
f.close()
# Save the matrices to files
f_trans = open('transition_matrices.dat','w')
f_rate  = open('rate_matrices.dat','w')
for i in range(N):
    for j in range(N):
        f_trans.write("({0:d},{1:d})  ".format(combined_transition_matrix[i,j][0], combined_transition_matrix[i,j][1]))
        f_rate.write( "({0:.4g},{1:.4g})  ".format(combined_rate_matrix[i,j][0]  , combined_rate_matrix[i,j][1]))
    f_trans.write("\n")
    f_rate.write("\n")
f_trans.close()
f_rate.close()

