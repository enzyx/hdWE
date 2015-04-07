from __future__ import print_function
import numpy
import constants
import analysis_operations

def reweightBinProbabilities(iterations, reweighting_range, workdir, jobname):
    '''
    Reweights the bin probabilities according to the steady state equations, using
    the mean rates over the last reweighting_range*len(iteration) iterations. 
    If the steady state equation matrix is singular, the bins added in the latest iterations are successively
    skipped for reweighting until the steady state equations can be solved.
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
    mean_rate_matrix= analysis_operations.getMeanRateMatrixWithConvergedOutrates(iterations, reweighting_range)
    #for j in range(len(mean_rate_matrix)):
    #    print(sum(mean_rate_matrix[j,:]))
    print('\nOld Bin Probabilities:', file=logfile)
    print(analysis_operations.meanBinProbabilities(iterations, iteration_counter, iteration_counter), file=logfile) 
    print('Using the last ' + str(iteration_range) + ' iterations for reweighting.', file=logfile)
    print('\nTotal Mean Rate Matrix:', file=logfile)
    print(mean_rate_matrix, file=logfile)
    for skip_number_of_latest_bins in range(0,iteration_range+1):
        last_bin_number=iterations[iteration_counter - skip_number_of_latest_bins].getNumberOfBins()
        skipped_bins = iterations[iteration_counter].getNumberOfBins() - last_bin_number

        print('Trying to solve steady state equations, skipping the bins added in the last ' + str(skip_number_of_latest_bins) + ' iteration(s) (skipping  ' + str(skipped_bins) + ' bin(s) ):', file=logfile)

        #delete entries for only in the last iteration newly created bins
        mean_rate_matrix=numpy.delete(mean_rate_matrix,numpy.s_[last_bin_number:iterations[iteration_counter].getNumberOfBins()],0)       
        mean_rate_matrix=numpy.delete(mean_rate_matrix,numpy.s_[last_bin_number:iterations[iteration_counter].getNumberOfBins()],1) 
        #get total Probability of the bins that will be reweighted
        prob_of_reweighted_bins = 0.0
        for i in range(0,last_bin_number):
            prob_of_reweighted_bins += iterations[-1].bins[i].getProbability()
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
                print('At least one Bin probablity would be zero.', file=logfile)
            else:
                #Assign the new probabilities to bins. Keep relative segment probabilities within bins.
                for i in range(0,len(bin_probs_from_rates)):
                    if iterations[-1].bins[i].getNumberOfSegments() > 0:
                        reweight_factor = 1.0 * bin_probs_from_rates[i] / iterations[-1].bins[i].getProbability()
                        for segment_loop in iterations[-1].bins[i]:
                            segment_loop.setProbability(segment_loop.getProbability() * reweight_factor)
                    else:
                        #TODO case not yet correctly handled
                        #pass
                        iterations[-1].bins[i].respawnSegmentFromReference(bin_probs_from_rates[i])
                print('Success.', file=logfile)            
                #print('Reduced Mean Rate Matrix:')
                #print(mean_rate_matrix) 
                print('Reweighted Bin Probabilities (Omitting latest ' + str(skipped_bins) + ' bins):', file=logfile)
                print(str(bin_probs_from_rates), file=logfile)
                return
        except:
            print('Singular rate matrix.', file=logfile)
 
    print('Singular rate matrix for skipping maximum number of latest bins. Skipping reweighting.', file=logfile)
    return
