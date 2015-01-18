import numpy

def meanRateMatrix(iterations, begin, end):
    """
    Calculates the mean values of the bin-to-bin rates over a given
    range of iterations  (end iteration is including) and returns the mean rate Matrix.
    Rates from or to bins that do not exist in a iteration are set to zero.
    It is assumed that bins are only added during iterating.
    """
    n_iterations = end - begin +1
    # Initialize mean rate matrix
    mean_rate_matrix = numpy.zeros([iterations[end].getNumberOfBins(), iterations[end].getNumberOfBins()], float)    
    
    # Sum rates
    for iterations_loop in iterations[begin:end + 1]:
        mean_rate_matrix_tmp = iterations_loop.RateMatrix()
        for i in range(0,len(mean_rate_matrix_tmp)):
            for j in range(0,len(mean_rate_matrix_tmp)):
                mean_rate_matrix[i,j] += mean_rate_matrix_tmp[i,j]
                
    # Divide by number of iterations            
    for i in range(0,len(mean_rate_matrix)):
        for j in range(0,len(mean_rate_matrix)): 
             mean_rate_matrix[i,j] =  mean_rate_matrix[i,j] / (n_iterations)
             
    return mean_rate_matrix
    
def meanBinProbabilities(iterations, begin, end):
    """
    Calculates the mean bin probabilities over a given
    range of iterations  (end iteration is including).
    It is assumed that bins are only added during iterating.
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
    rate_matrix = numpy.transpose(rate_matrix)
    for i in range(0,len(rate_matrix)):
        total_flux_out    = sum(rate_matrix[:,i])
        rate_matrix[i,i]  +=  - total_flux_out
    
    # Make the first equation 1.0 * bin[0].probability = 1.0 to set the first probability to a nonzero value.
    # TODO first matrix row in a way that it includes normalization?
    rate_matrix[0,0] = 1.0
    for i in range(1,len(rate_matrix)):
        rate_matrix[0,i] = 0.0
    # Use numpy to solve the system of linear equations 
    bin_probabilities = numpy.linalg.solve(rate_matrix, zeros_vector)
        
    # Norm
    total_probability = sum(bin_probabilities)
    for i in range(0,len(bin_probabilities)):
        bin_probabilities[i] = bin_probabilities[i] / total_probability
        
    return bin_probabilities