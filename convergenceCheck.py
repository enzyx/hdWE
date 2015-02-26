from __future__ import print_function
import numpy as np
import constants
import analysis_operations

def checkOutratesForConvergence(iterations, current_iteration, convergence_range, max_rate_rmsf):
    '''
        checks outgoing rates of all bins for convergence. 
        If all outgoing rates are converged the bin is
        flagged as having converged outrates.
        
        Convergence is achieved if the rmsf of the
        last <convergence_range> iterations is below threshold.
    '''
    iteration_counter = len(iterations) - 1
    if iteration_counter > convergence_range:
        rate_matrices = []
        """
           calculate mean rate matrices averaged 
           over the last <reweighting_range>:=rr, rr-1, rr-2 etc
           iterations
        """
        # number of bins that exist long enough to have converged rates
        n_bins_to_check = iterations[iteration_counter - \
        							  convergence_range].getNumberOfBins()
        #~ print ("checking {} bins".format(n_bins_to_check))

        # calculate mean rate matrices
        for i in range(convergence_range):
            rate_matrices.append(iterations[iteration_counter - convergence_range + i].RateMatrix())
            #~ print(rate_matrices[-1])
        	
        #check for rate-converged bins        
        for bin_index in range(n_bins_to_check): # for every bin that exists long enough to have converged rates
            if not current_iteration.bins[bin_index].isConverged():
                b_all_rates_exist = True
                # get all rates
                # rates is an array of rates to all bins
                all_bin_rates = []
                # get rates to all bins that in the last matrix.
                for target_bin_index in range(len(rate_matrices[-1][bin_index])):
                    all_bin_rates.append(np.array([ m[bin_index][target_bin_index] \
                                                    for m in rate_matrices\
                                                    if len(m[bin_index]) > target_bin_index ]))
                    if  np.amax(all_bin_rates[-1]) > constants.num_boundary and len(all_bin_rates[-1]) != convergence_range:                                                     
                        # some non-zero rates haven't been measured long enough
                        b_all_rates_exist = False
                        
                # debug output
                #~ if bin_index == 0:
                    #~ for rates in all_bin_rates:
                        #~ mean = np.mean(rates)
                        #~ rmsf = np.sqrt(np.mean(np.square(rates - mean)))
                        #~ if mean < constants.num_boundary:
                            #~ if rmsf < constants.num_boundary:
                                #~ relative_rmsf = 0.0
                            #~ else:
                                #~ #TODO: what to do with small mean but larger rmsf? is that even possible?
                                #~ relative_rmsf = rmsf
                        #~ else:
                            #~ relative_rmsf = rmsf/mean
                        #~ print (rates, ", mean", mean, ", rmsf", rmsf, ", relative_rmsf", relative_rmsf)
                
                if b_all_rates_exist:
                    #~ print ("all rates present. checking convergence of rates for bin", bin_index)
                    b_converged = True    
                    for target_bin_index,rates in enumerate(all_bin_rates):
                        mean = np.mean(rates)
                        rmsf = np.sqrt(np.mean(np.square(rates - mean)))
                        if mean < constants.num_boundary:
                            if rmsf < constants.num_boundary:
                                relative_rmsf = 0.0
                            else:
                                # in the exotic case that the mean is basically zero 
                                # but the rmsf isn't we don't consider it converged
                                relative_rmsf = max_rate_rmsf
                                print ("relative_rmsf({}) > num_boundary but mean{} < num_boundary!.".format(rmsf, mean) + 
                                              "Case not handled correctly yet.")
                        else:
                            relative_rmsf = rmsf/mean
                        #~ print ("rmsf/mean for {bin}->{binout}: {rrmsf}".format(bin=bin_index, binout=target_bin_index, rrmsf=relative_rmsf))

                        if relative_rmsf >= max_rate_rmsf:
                            b_converged = False
                            break

                    if b_converged:
                        current_iteration.bins[bin_index].setConverged(True)
                        print ("\nBin {} is converged".format(bin_index))
                					
        #~ print ("rates[0]:")				
        #~ print (rates[0])
    return
