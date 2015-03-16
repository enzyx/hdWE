import numpy
from copy import deepcopy
import constants
    
def meanRateMatrix(iterations, begin, end):
    """
    Calculates the mean values of the bin-to-bin rates over a given
    range of iterations  (end iteration is including) and returns the mean rate Matrix.
    The mean value of outrates of bins that do not exists over the full iteration range is calculated
    only over the existing iterations in order to not underestimate them.
    """
    # Initialize mean rate matrix
    mean_rate_matrix = numpy.zeros([iterations[end].getNumberOfBins(), iterations[end].getNumberOfBins()], float)  
    bin_exists_since = numpy.zeros([iterations[end].getNumberOfBins()], int)  
    
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
    
def getMeanRateMatrixWithConvergedOutrates(iterations, reweighting_range):
    """
    Calculates the mean Rate Matrix over a range of iterations going back from
    the iteration at which the outrate was converged.
    The iteration_range is adapted to the iteration the bins was converged.
    Thus, once a bin is converged, the outrates stay the same.
    (In this way the outrates of converged bins can be recalculated from the information
    in the logfile and they dont have to be explicitly saved.)
    """
    end   = len(iterations - 1)
    begin = int( end * (1.0 - reweighting_range) )
    mean_rate_matrix = meanRateMatrix(iterations, begin, end)
    # get iterations of convergence
    iteration_of_convergence = numpy.zeros([iterations[end].getNumberOfBins()],int)
    for iteration_loop in iterations:
        for bin_loop in iteration_loop:
            if not bin_loop.isConverged():
                iteration_of_convergence[bin_loop.getId()] = iteration_loop.getId()
    for bin_loop in iterations[-1]:
        if bin_loop.isConverged():
            end_tmp   = iteration_of_convergence[bin_loop.getId()]
            begin_tmp = int( end * (1.0 - reweighting_range) )
            temp_mean_rate_matrix = meanRateMatrix(iterations,
                                    begin_tmp,
                                    end_tmp)
            for i in range(len(temp_mean_rate_matrix[0,:])):
                mean_rate_matrix[bin_loop.getId(),i] = temp_mean_rate_matrix[bin_loop.getId(),i]
            for i in range(len(temp_mean_rate_matrix[0,:]),len(mean_rate_matrix)):
                mean_rate_matrix[bin_loop.getId(),i] = 0.0
           
    return mean_rate_matrix
        
        
    
