"""
Calculates the total Rates in and out of a bin.
"""
from __future__ import print_function
import numpy
import constants
from logger import Logger
import math
from amber_module import MD_module
import argparse  
import analysis_operations

###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-c', '--conf', dest="input_md_conf", 
                    type=str, default = 'hdWE.conf',
                    help="MD-Software configuration file")
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0,
                    help="First iteration to use for PMF calculation.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Last iteration to to use for PMF calculation.")  
parser.add_argument('-i', '--bin_index', dest="bin_index",
                    type=int, default=0,
                    help="Bin index.") 
                    
# Initialize
print('\033[1mAnalyzing Rates.\033[0m')      
args = parser.parse_args()
md_module = MD_module(args.input_md_conf,debug=False)

#get the actual Iteration from logger module
logger = Logger(args.logfile, append = True)
iterations = logger.loadIterations(0, args.last_iteration, bCheckFiles = False)
logger.close()

first_it = args.first_iteration
last_it  = iterations[-1].getId()
bin_index = args.bin_index

rate_matrix = analysis_operations.meanRateMatrix(iterations, first_it, last_it)

outrate = numpy.sum(rate_matrix[bin_index,:])
inrate  = numpy.sum(rate_matrix[:,bin_index])

#subtract rate from bin to itself
outrate -= rate_matrix[bin_index,bin_index]
inrate  -= rate_matrix[bin_index,bin_index]

print('Total Rate Out:' + str(outrate))
print('Total Rate In: ' + str(inrate))
print('-> kT*ln(outrate/inrate):' + str(constants.kT * math.log(outrate / inrate)))



