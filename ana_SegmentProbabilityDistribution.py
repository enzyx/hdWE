#!/usr/bin/env python2
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
Analyzes the distribution of probabilities of the segments of a bin.
"""
from __future__ import print_function
import numpy
import constants
import sys
import segment
from logger import Logger
from math import log
from amber_module import MD_module
import argparse  

###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-c', '--conf', dest="input_md_conf", nargs='?',
                    default=False, 
                    help="Optional MD-Software configuration file")
parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    default="log", metavar="DIR", required=True,
                    help="The log directory")
parser.add_argument('-i', '--it', dest="iteration",
                    type=int, default=0,
                    help="Iteration to analyize.")       
parser.add_argument('-o', '--output', dest="output_path", 
                    type=str, default='ana_SegmentProbabilityDistribution.output',
                    help="Output filename")  
                    
# Initialize
print('\033[1mAnalyzing segment probabilities\033[0m')      
args = parser.parse_args()

#get the actual Iterations from logger module
logger = Logger(args.logdir)
iterations = logger.loadIterations(begin=args.iteration, end=args.iteration)
iteration  = iterations[-1]

statistics        = numpy.zeros([iteration.getNumberOfBins(),5], float)
probabilities_tmp = numpy.zeros([iteration.bins[0].getNumberOfSegments()]) 
for i in range(iteration.getNumberOfBins()):
    if iteration.bins[i].getNumberOfSegments()>0:
        for j in range(iteration.bins[i].getNumberOfSegments()):
            probabilities_tmp[j] = iteration.bins[i].segments[j].probability        
        statistics[i,0] = max(probabilities_tmp)  / min(probabilities_tmp) 
        statistics[i,1] = max(probabilities_tmp) 
        statistics[i,2] = min(probabilities_tmp)
        statistics[i,3] = numpy.mean(probabilities_tmp)
        statistics[i,4] = numpy.std(probabilities_tmp)
    
header = ('seg. probabilities: max_prob/min_prob max_prob min_prob mean std')
numpy.savetxt(args.output_path, statistics, header=header)



