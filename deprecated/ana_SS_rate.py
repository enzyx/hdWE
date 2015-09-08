#!/usr/bin/env python2
"""
Calculates the rate from a hdWE steady state run.
"""
from __future__ import print_function
import numpy
import scipy.stats
from lib.logger import Logger
from lib.functions_ana_general import autocorrelation_function
import argparse 
#import scikits.bootstrap
 

###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-l', '--log', dest="logdir",
                    required=True, 
                    help="hdWE log directory")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0,
                    help="First iteration to use.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Last iteration to to use.")  
parser.add_argument('-o', '--output', dest="output_path", 
                    type=str, default='ana_SS_rate',
                    help="Output filename")  
               
            
               
# Initialize
print('\033[1mCalculating Rate.\033[0m')      
args = parser.parse_args()
logger = Logger(args.logdir)
iterations = logger.loadIterations(begin=args.first_iteration, end=args.last_iteration)

# get the recycled probability and number of recycled segments in each iteration
# recycled segments are those that reached end state bins

# Mean Recycled Probability
recycled_probability = []
N_recycled_segments  = 0
for this_iteration in iterations:
    if len(this_iteration.getEndStateBins()) > 0:
        for this_bin in this_iteration.getEndStateBins():
                recycled_probability.append(this_bin.getInitialProbability())
                N_recycled_segments += this_bin.getNumberOfInitialSegments()
    else:
        recycled_probability.append(0.0)
# Get probability in start state
p_start = 0.0
for this_iteration in iterations:
    if len(this_iteration.getStartStateBins()) > 0:
        for this_bin in this_iteration.getStartStateBins():
            p_start += this_bin.getProbability()
p_start /= len(iterations)
print("Mean start state probability = {:5.4e}".format(p_start))

meanCI_rate = scipy.stats.bayes_mvs(recycled_probability, 0.95)[0]

# Cumulative Mean Recycled Probability
recycled_probability_cumulative = numpy.zeros([len(recycled_probability),3])
for i in range(2,len(recycled_probability)):
    recycled_probability_cumulative_tmp = scipy.stats.bayes_mvs(recycled_probability[0:i], 0.95)[0]
    recycled_probability_cumulative[i,0] = recycled_probability_cumulative_tmp[0]
    recycled_probability_cumulative[i,1] = recycled_probability_cumulative_tmp[1][0]    
    recycled_probability_cumulative[i,2] = recycled_probability_cumulative_tmp[1][1]


# Autocorrelation Function
autocorr = autocorrelation_function(recycled_probability)

# print
numpy.savetxt(args.output_path+'.rate', recycled_probability)
numpy.savetxt(args.output_path+'.rate_cumulative', recycled_probability_cumulative) 
numpy.savetxt(args.output_path+'.autocorrelation', autocorr)

print(' Number of recycled segments: {:9d}'.format( N_recycled_segments ) )
print(' Rate: {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
      mean = meanCI_rate[0],
      lCI  = meanCI_rate[1][0],
      uCI  = meanCI_rate[1][1]))

