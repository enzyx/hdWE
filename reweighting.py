from __future__ import print_function
import numpy
import constants
import analysis_operations

def reweightBinProbabilities(iterations, iteration_range):
    '''
    Reweights the bin probabilities according to the steady state equations, using
    the mean rates over the last iteration_range iterations. 
    If a new bin was created within the last iteration_range iterations, this bin is omitted for the reweighting,
    because no reliable out-rate is yet available. 
    The total probablity of the reweighted bins is not changed.
    The bin probablities of skipped bins are not changed.
    '''
    iteration_counter = len(iterations) - 1
    if iteration_counter >= iteration_range:
        print('\nOld Bin Probabilities:')
        print(analysis_operations.meanBinProbabilities(iterations, iteration_counter, iteration_counter))         
        #get Rate Matrix, averaged over the last iteration_range iterations
        mean_rate_matrix= analysis_operations.getMeanRateMatrixWithConvergedOutrates(iterations, iteration_counter - iteration_range, iteration_counter)
        print('\nTotal Mean Rate Matrix:')
        print(mean_rate_matrix)         
        #delete entries for only in the last iteration newly created bins
        last_bin_number=iterations[iteration_counter - iteration_range].getNumberOfBins()
        skipped_bins = iterations[iteration_counter].getNumberOfBins() - last_bin_number
        mean_rate_matrix=numpy.delete(mean_rate_matrix,numpy.s_[last_bin_number:iterations[iteration_counter].getNumberOfBins()],0)       
        mean_rate_matrix=numpy.delete(mean_rate_matrix,numpy.s_[last_bin_number:iterations[iteration_counter].getNumberOfBins()],1) 
        print('\nReduced Mean Rate Matrix:')
        print(mean_rate_matrix)         
        #get total Probability of the bins that will be reweighted
        prob_of_reweighted_bins = 0.0
        for i in range(0,last_bin_number):
            prob_of_reweighted_bins += iterations[iteration_counter].bins[i].getProbability()
        try:
            #Try to solve the steady state equations:
            bin_probs_from_rates = analysis_operations.BinProbabilitiesFromRates(mean_rate_matrix)
            #Rescale the normalized bin probabilities from the state state equations
            #to the total probablity they had before:
            for i in range(0,len(bin_probs_from_rates)):
                bin_probs_from_rates[i] = bin_probs_from_rates[i] * prob_of_reweighted_bins
            #Skip reweighting if a bin would become 0
            #(this could lead to trajectories with 0 probablity assigned)
            if numpy.min(bin_probs_from_rates) <= constants.num_boundary:
                print('At least one new probability would be 0. Skipping reweighting.')
                return
            #Assign the new probabilities to bins. Keep relative segment probabilities within bins.
            for i in range(0,len(bin_probs_from_rates)):
                reweight_factor = bin_probs_from_rates[i] / iterations[iteration_counter].bins[i].getProbability()
                if iterations[iteration_counter].bins[i].getNumberOfSegments() > 0:
                    for segment_loop in iterations[iteration_counter].bins[i]:
                        segment_loop.setProbability(segment_loop.getProbability() * reweight_factor)
                else:
                    iterations[iteration_counter].bins[i].respawnSegmentFromReference(bin_probs_from_rates[i])
                        
            if skipped_bins==0:
                print('Reweighted Bin Probabilities:')                
            else:
                print('Reweighted Bin Probabilities (Omitting latest ' + str(skipped_bins) + ' bins):')
            print(bin_probs_from_rates)            
        
        except:
            print('Singular rate matrix. Skipping reweighting.')
 
    #else:
    #    print('Iteration progress smaller then reweighting iteration range. Skipping reweighting.')

    return