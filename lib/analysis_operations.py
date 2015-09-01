import numpy
from copy import deepcopy
import lib.constants as constants
from lib.segment import Segment
    
def meanRateMatrix(iterations, begin=0, end=-1):
    """
    Calculates the mean values of the bin-to-bin rates over a given
    range of iterations  (end iteration is including) and returns the mean rate Matrix.
    The mean value of outrates of bins that do not exists over the full iteration range is calculated
    only over the existing iterations in order to not underestimate them.
    """
    # Accept default values for begin and end
    if end == -1:
        end = len(iterations) - 1 
    
    # Initialize mean rate matrix
    N = iterations[end].getNumberOfBins()
    mean_rate_matrix = numpy.zeros([N,N], float)  
    bin_exists_since = numpy.zeros([N], int)  
    
    # Sum rates
    for iterations_loop in iterations[begin:end + 1]:
        mean_rate_matrix_tmp = iterations_loop.RateMatrix()
        for i in range(0,len(mean_rate_matrix_tmp)):
            if max(mean_rate_matrix_tmp[i,:]) > constants.num_boundary:
                bin_exists_since[i] += 1
                for j in range(0,len(mean_rate_matrix_tmp)):
                    mean_rate_matrix[i,j] += mean_rate_matrix_tmp[i,j]
                
    # Divide by number of iterations            
    for i in range(0,len(mean_rate_matrix)):
        if bin_exists_since[i] > 0:
            for j in range(0,len(mean_rate_matrix)):
                mean_rate_matrix[i,j] =  mean_rate_matrix[i,j] / bin_exists_since[i]
             
    return mean_rate_matrix

def cumulativeTransitionMatrix(iterations, begin=0, end=-1):
    """
    Return the cumulative Transition matrix of iterations[begin, end] 
    (end is inculded).
    """
    # Accept default values for begin and end
    if end == -1:
        end = len(iterations) - 1 
    
    N = iterations[end].getNumberOfBins()
    cumulative_transition_matrix = numpy.zeros((N,N), int)  
    for this_iteration in iterations[begin:end + 1]:
        this_transition_matrix = this_iteration.TransitionMatrix()
        M = len(this_transition_matrix)
        for i in range(M):
            for j in range(M):
                cumulative_transition_matrix[i,j] += this_transition_matrix[i,j]
    return cumulative_transition_matrix
    
def meanBinProbabilities(iterations, begin, end):
    """
    Calculates the mean bin probabilities over a given
    range of iterations  (end iteration is including).
    """
    n_iterations = end - begin +1
    # Initialize bin probabilities
    bin_probabilities = numpy.zeros([iterations[end].getNumberOfBins()], float)    
    
    # Sum probabilities
    for iterations_loop in iterations[begin:end + 1]:
        for i in range(0,len(iterations_loop.bins)):
            bin_probabilities[i] += iterations_loop.bins[i].getProbability()
                
    # Divide by number of iterations            
    for i in range(0,len(bin_probabilities)):
        bin_probabilities[i] = bin_probabilities[i] / n_iterations
             
    return bin_probabilities
    
def BinProbabilitiesFromRates(rate_matrix):
    """
    Solve the steady state condition (bin_flux_in = bin_flux_out) for a given rate matrix
    and return the corresponding bin probabilities.
    """
    # Initialize matrix and vectors
    zeros_vector    = numpy.zeros([len(rate_matrix)], float)
    zeros_vector[0] = 1.0 # condition that the probability of the first bin is 1.
    # Matrix corresponding to the system of linear equations: 
    # Transpose rate matrix and subtract total flux out of bin from the diagonal entries
    equation_matrix = deepcopy(numpy.transpose(rate_matrix))
   
    for i in range(0,len(equation_matrix)):
        total_flux_out        = sum(equation_matrix[:,i])
        equation_matrix[i,i]  +=  - total_flux_out
    # Make the first equation 1.0 * bin[0].probability = 1.0 to set the first probability to a nonzero value.
    # TODO first matrix row in a way that it includes normalization?
    equation_matrix[0,0] = 1.0
    for i in range(1,len(equation_matrix)):
        equation_matrix[0,i] = 0.0
    # Use numpy to solve the system of linear equations 
    bin_probabilities = numpy.linalg.solve(equation_matrix, zeros_vector)
        
    # Norm
    total_probability = sum(bin_probabilities)
    for i in range(0,len(bin_probabilities)):
        bin_probabilities[i] = bin_probabilities[i] / total_probability
    
    return bin_probabilities

