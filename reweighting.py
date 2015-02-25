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

def checkOutratesForConvergence(iterations, convergence_range, max_rate_rmsf):
    '''
        checks outgoing rates of all bins for convergence. 
        If all outgoing rates are converged the bin is
        flagged as having converged outrates.
        
        Convergence is achieved if the rmsf of the
        last <convergence_range> iterations is below threshold.
    '''
    iteration_counter = len(iterations)
    if iteration_counter > convergence_range:
        rate_matrices = []
        """
           calculate mean rate matrices averaged 
           over the last <reweighting_range>:=rr, rr-1, rr-2 etc
           iterations
        """
        #~ numpy.set_printoptions(precision=3)
        n_bins_to_check = iterations[iteration_counter - \
        							  convergence_range].getNumberOfBins()
        #~ print ("max bins mutable = {}".format(n_bins_to_check))

        # calculate mean rate matrices
        #~ print("\nrate matrices:")
        for i in range(convergence_range):
            rate_matrices.append(analysis_operations.meanRateMatrix(iterations,
        											 iteration_counter - i - 1,
        											 iteration_counter - i - 1))
            #~ print(rate_matrices[-1])
        	
        #check for rate-converged bins
        rates=[[numpy.zeros(convergence_range)\
            for j in range(n_bins_to_check)]\
             for i in range(n_bins_to_check)]
        for bin_index in range(n_bins_to_check): # for every bin
            for matrix_index in range(2): # for every matrix
                for target_bin_index in range(n_bins_to_check):
                    if target_bin_index != bin_index:
                        rates[bin_index][target_bin_index][matrix_index] = \
                            rate_matrices[matrix_index][bin_index][target_bin_index]
                            
                

            # check convergency
            if len(rates[bin_index]) > 0:
                bConverged = True
                for target_bin_index in range(n_bins_to_check):
                    mean = numpy.mean(rates[bin_index][target_bin_index])
                    rmsf = numpy.sqrt(numpy.mean(numpy.square(\
                            rates[bin_index][target_bin_index] - mean)))
                    #~ print ("rmsf bin {bin} to {binout}: {rmsf}".format(bin=bin_index, binout=target_bin_index, rmsf=rmsf))
                    if rmsf >= max_rate_rmsf:
                        bConverged = False
                        break
                if bConverged:
                    iterations[-1].bins[bin_index].set_outrates_converged(bConverged)
                    print ("\n Bin {} is converged".format(bin_index))
                					
        #~ print ("rates[0]:")				
        #~ print (rates[0])
    return
