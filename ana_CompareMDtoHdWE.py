#!/usr/bin/env python2
"""
Compare MD and hdWE rama bin properties 
"""
from __future__ import print_function
import numpy
from logger import Logger
import argparse  


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
parser.add_argument('-m', '--md-probabilities', type=str, dest="md_probs", 
                    metavar="FILE", default=False, nargs="?", const="md_probabilities.dat",
                    help="Write bin data sorted by MD probabilities.")
parser.add_argument('-w', '--hdWE-probabilities', type=str, dest="hdWE_probs", 
                    metavar="FILE", default=False, nargs="?", const="hdWE_probabilities.dat",
                    help="Write bin data sorted by hdWE probabilities.")

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

md_file = open(args.md_rama_bins,'r')
nframes = 0
for line in md_file.readlines():
    nframes += 1.0
    frame = line.strip().split()
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

# Sort by p_md and write to output file if requested by user
if args.md_probs:
    f = open(args.md_probs, 'w')
    f.write("# rama_id   p_md   p_hdWE   bin_id\n")
    for b in sorted(bin_data, key=lambda e:e['p_md'], reverse=True):
        f.write("{0:s} {1:15e} {2:15e} {3:d}\n".format(b["rama_id"], b["p_md"], b["p_hdWE"], b["id"]))
    f.close()

# Sort by p_hdWE and write to output file if requested by user
if args.hdWE_probs:
    f = open(args.hdWE_probs, 'w')
    f.write("# rama_id   p_md   p_hdWE   bin_id")
    for b in sorted(bin_data, key=lambda e:e['p_hdWE'], reverse=True):
        f.write("{0:s} {1:15e} {2:15e} {3:d}\n".format(b["rama_id"], b["p_md"], b["p_hdWE"], b["id"]))
    f.close()