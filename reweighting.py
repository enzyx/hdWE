from __future__ import print_function
import numpy
import analysis_operations

def reweightBinProbabilities(iterations, iteration_range):
    '''
    Reweights the bin probabilities according to the steady state equations, using
    the mean rates over the last iteration_range iterations. 
    If a new bin was created within the last iteration_range iterations, this bin is omitted for the reweighting,
    because no reliable out-rate is yet available. 
    The reweighted and not reweighted bin probabilities are then normalized.
    '''
    iteration_counter = len(iterations) - 1
    if iteration_counter >= iteration_range:
        print('\nOld Bin Probabilities:')
        print(analysis_operations.meanBinProbabilities(iterations, iteration_counter, iteration_counter))         
       #get Rate Matrix, averaged over the last iteration_range iterations
        mean_rate_matrix= analysis_operations.meanRateMatrix(iterations, iteration_counter - iteration_range, iteration_counter)
        #delete entries for only in the last iteration newly created bins
        last_bin_number=iterations[iteration_counter - iteration_range].getNumberOfBins()
        mean_rate_matrix=numpy.delete(mean_rate_matrix,numpy.s_[last_bin_number:iterations[iteration_counter].getNumberOfBins()],0)       
        mean_rate_matrix=numpy.delete(mean_rate_matrix,numpy.s_[last_bin_number:iterations[iteration_counter].getNumberOfBins()],1) 
        try:
            #Try to solve the steady state equations:
            bin_probs_from_rates = analysis_operations.BinProbabilitiesFromRates(mean_rate_matrix)
            #add old non-reweighted bin probabilities
            new_bin_probs = numpy.zeros([iterations[iteration_counter].getNumberOfBins()],float)
            skipped_bins = len(new_bin_probs) - last_bin_number 
            new_bin_probs[0:last_bin_number]=bin_probs_from_rates
            for i in range(last_bin_number,len(new_bin_probs)):
                new_bin_probs[i] = iterations[iteration_counter].bins[i].getProbability()
            #Normalize the complete array of new bin probabilities:
            total_prob_tmp = numpy.sum(new_bin_probs)
            for i in range(0,len(new_bin_probs)):
                new_bin_probs[i] = new_bin_probs[i] / total_prob_tmp
            #Assign the new probabilities to bins. Keep relative segment probabilities within bins.
            for bin_loop in iterations[iteration_counter].bins:
                for segment_loop in bin_loop:
                    reweight_factor = new_bin_probs[bin_loop.getId()] / bin_loop.getProbability()
                    if bin_loop.getNumberOfSegments()>0.0:
                         segment_loop.probability = segment_loop.probability * reweight_factor
            if skipped_bins==0:
                print('Reweighted Bin Probabilities:')                
            else:
                print('Reweighted Bin Probabilities (approximated: skipped latest ' + str(skipped_bins) + ' bins):')
            print(new_bin_probs)            
            
        except:
            print('Singular rate matrix. Skipping reweighting.')
 
    #else:
    #    print('Iteration progress smaller then reweighting iteration range. Skipping reweighting.')

    return