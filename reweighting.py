from __future__ import print_function
import numpy
import constants
import analysis_operations

def reweightBinProbabilities(iterations, reweighting_range, workdir, jobname):
    '''
    Reweights the bin probabilities according to the steady state equations, using
    the mean rates over the last reweighting_range*len(iteration) iterations. 
    Bins with no outrates to other bins are excluded from the reweighting.
    The total probablity of the reweighted bins is not changed.
    The bin probablities of skipped bins are not changed.
    '''
    
    iteration_range = int(len(iterations) * reweighting_range)
    iteration_counter = len(iterations) - 1
    logfile=open(workdir+jobname+'-log/reweighting'+str(iteration_counter).zfill(5),'w+')
    #check for empty bins    
    #for bin_loop in iterations[-1].bins:
    #    if bin_loop.getNumberOfSegments() < 1:
    #        print('At least one bin is empty. Skipping reweighting.', file=logfile)
    #        return
    #get Rate Matrix, averaged over the last iteration_range iterations
    mean_rate_matrix       = analysis_operations.getMeanRateMatrixWithConvergedOutrates(iterations, reweighting_range)
    bin_probabilities = iterations[-1].getBinProbabilities()
    #for j in range(len(mean_rate_matrix)):
    #    print(sum(mean_rate_matrix[j,:]))
    print('\nOld bin probabilities:', file=logfile)
    print(bin_probabilities, file=logfile)
    print('Using the last ' + str(iteration_range) + ' iterations for reweighting.', file=logfile)
    print('\nTotal Mean Rate Matrix:', file=logfile)
    print(mean_rate_matrix, file=logfile)
    print('Old chi square (equilibriums condition): {0:g}'.format(
          calcChiSquare(bin_probabilities, mean_rate_matrix)), file=logfile)
    #check for bins with no outrates or only outrate of 1 to itself
    keep_bin_index   = numpy.zeros([0], int)
    delete_bin_index = numpy.zeros([0], int)
    for i in range(0,len(mean_rate_matrix)):
        if max(mean_rate_matrix[i,:]) > 0.0 and max(mean_rate_matrix[i,:]) < 1.0:
            keep_bin_index = numpy.append(keep_bin_index, i)
        else: 
            delete_bin_index = numpy.append(delete_bin_index, i)  
    #reduce mean_rate_matrix
    reduced_mean_rate_matrix = numpy.delete(mean_rate_matrix, delete_bin_index, 0)
    reduced_mean_rate_matrix = numpy.delete(reduced_mean_rate_matrix, delete_bin_index, 1)
    
    #get total Probability of the bins that will be reweighted
    prob_of_reweighted_bins = 0.0
    for i in range(0,len(keep_bin_index)):
        prob_of_reweighted_bins += iterations[-1].bins[keep_bin_index[i]].getProbability()
        
    print('Trying to solve steady state equations for the reduced rate matrix:', file=logfile)
    #Try to solve the steady state equations:
    try:
        bin_probs_from_rates = analysis_operations.BinProbabilitiesFromRates(reduced_mean_rate_matrix)
        #Rescale the normalized bin probabilities from the steady state equations
        #to the total probability they had before:
        for i in range(0,len(keep_bin_index)):
            bin_probs_from_rates[i] = bin_probs_from_rates[i] * prob_of_reweighted_bins

        #Skip reweighting if a bin would become 0
        #(this could lead to trajectories with 0 probability assigned)
        if numpy.min(bin_probs_from_rates) <= constants.num_boundary:
            print('At least one Bin probability would be zero.', file=logfile)
            print('At least one Bin probablity would be zero.')
        else:
            #Assign the new probabilities to bins. Keep relative segment probabilities within bins.
            for i in range(0,len(keep_bin_index)):
                if iterations[-1].bins[keep_bin_index[i]].getNumberOfSegments() > 0:
                    reweight_factor = 1.0 * bin_probs_from_rates[i] / iterations[-1].bins[keep_bin_index[i]].getProbability()
                    for segment_loop in iterations[-1].bins[keep_bin_index[i]]:
                        segment_loop.setProbability(segment_loop.getProbability() * reweight_factor)
                else:
                    print('Created new segments in bin ' +str(keep_bin_index[i])+ ' by reweighting', file = logfile)            
                    iterations[-1].bins[keep_bin_index[i]].respawnSegmentFromReference(bin_probs_from_rates[i])
            print('Success.', file=logfile)    
            print('Skipped bins ' + str(delete_bin_index), file=logfile)
            print('Reweighted Bin Probabilities:', file=logfile)
            print(str(bin_probs_from_rates), file=logfile)
            
            print('     Reweighting: Success. Skipped bins: ' + str(delete_bin_index))
            bin_probabilities = iterations[-1].getBinProbabilities()
            print('New chi square (equilibriums condition): {0:g}'.format(
                  calcChiSquare(bin_probabilities, mean_rate_matrix)), file=logfile)
            return
    except:
        print('Singular rate matrix.', file=logfile)
        print('\x1b[31m     Singular rate matrix.\x1b[0m')
    return

def calcChiSquare(P, K):
    """
    Calculates the chi squared deviation from the equilibrium condition 
    There are N*(N-1)/2 equilibriums conditions P_i*k_ji = P_j*k_ij
    @param P: Vector of bin probabilities
    @param K: Rate matrix k_ij 
    """ 
    chi2 = 0.0
    N = len(P)
    for i in range(0,N-1):
        for j in range(i+1,N):
            chi2 += (P[i] * K[i][j] - P[j]*K[j][i])**2
    return chi2/(N)**2
