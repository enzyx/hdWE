#!/usr/bin/python2
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
from __future__ import print_function
from logger import Logger
import analysis_operations
import argparse 
import numpy
import matplotlib.pyplot as plt

parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    default="hdWE-log", metavar="DIR",
                    help="The logfile for reading and writing")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0,
                    help="First iteration to use for calculation.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Last iteration to to use for calculation.") 
parser.add_argument('-o', '--output', dest="output_filename",
                    type=str, default='Equilibrium',
                    help="Output Filename.")                    
                    
# Initialize
print('\033[1mAnalyzing Bins and Rates.\033[0m')      
args = parser.parse_args()


#get the actual Iteration from logger module
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.first_iteration, args.last_iteration)
n_iterations = args.last_iteration - args.first_iteration + 1 

n_bins = iterations[-1].getNumberOfBins()
diff_matrix = numpy.zeros([n_bins, n_bins],float)
for i in range(0,len(diff_matrix[:,0])):
    for j in range(0,len(diff_matrix[0,:])):
        diff_matrix[i,j] = 'inf'

rates           = analysis_operations.meanRateMatrix(iterations, 0, len(iterations)-1 ) 
probabilities   = analysis_operations.meanBinProbabilities(iterations, 0, len(iterations)-1 )


for i in range(0,n_bins):
    for j in range(0,n_bins):
        if i-j<0:
            if rates[j,i] > 0.001  and rates[i,j] > 0.001 and probabilities[i] > 0.0 and probabilities[j] > 0.0:
                diff_matrix[i,j] = (rates[i,j]*probabilities[i]) / (rates[j,i]*probabilities[j])
                if diff_matrix[i,j] < 1.0:
                    diff_matrix[i,j] = 1.0 / diff_matrix[i,j]
                if diff_matrix[i,j] > 10.0:
                    diff_matrix[i,j] = 10.0
        
#Plot
fig=plt.figure(figsize=(5,5))
plt.title('Equilibrium Condition')
plt.xlabel('# bin')
plt.ylabel('# bin')
plt.imshow(diff_matrix, interpolation='none',origin='lower')
cbar=plt.colorbar()
plt.jet()
cbar.set_label('max(pifij/pjfji,pjfji/pifij)')
plt.savefig(args.output_filename+'.pdf',format='pdf',dpi=300) 


