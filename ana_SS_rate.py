#!/usr/bin/env python2
"""
Calculates the rate from a hdWE steady state run.
"""
from __future__ import print_function
import numpy
import scipy.stats
from logger import Logger
import argparse 
import scikits.bootstrap
 

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
                    type=str, default='ana_SS_rate.output',
                    help="Output filename")  
               
def autocorr(x, t):
    return numpy.corrcoef(numpy.array([x[0:len(x)-t], x[t:len(x)]]))               
               
# Initialize
print('\033[1mCalculating Rate.\033[0m')      
args = parser.parse_args()

# get the Iterations from logger module
logger = Logger(args.logdir)
iterations = logger.loadIterations(begin=args.first_iteration, end=args.last_iteration)

# get the recycled probability in each iteration, which is the probability
# that reached end state bins
recycled_probability = []
N_recycled_segments  = 0

for this_iteration in iterations:
    if len(this_iteration.getEndStateBins()) > 0:
        for this_bin in this_iteration.getEndStateBins():
                recycled_probability.append(this_bin.getInitialProbability())
                N_recycled_segments += this_bin.getNumberOfInitialSegments()
    else:
        recycled_probability.append(0.0)

auto = []
for i in range(len(recycled_probability)/2):
    auto.append(autocorr(recycled_probability, i)[0,1])

numpy.savetxt('auto.dat', auto)

meanCI_rate = scipy.stats.bayes_mvs(recycled_probability, 0.95)[0]
#print( scikits.bootstrap.ci(recycled_probability, numpy.mean, alpha = 0.05) )

recycled_probability_cumulative = []
for i in range(1,len(recycled_probability)):
    recycled_probability_cumulative.append(numpy.mean(recycled_probability[0:i]))


# print
numpy.savetxt(args.output_path, recycled_probability)
numpy.savetxt(args.output_path+'_cumulative', recycled_probability_cumulative)

print(' Number of recycled segments: {:9d}'.format( N_recycled_segments ) )
print(' Rate: {mean:5.3e}   CI: {lCI:5.3e} to {uCI:5.3e}'.format(
      mean = meanCI_rate[0],
      lCI  = meanCI_rate[1][0],
      uCI  = meanCI_rate[1][1]))

