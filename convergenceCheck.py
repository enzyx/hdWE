from __future__ import print_function
import sys
import numpy as np
import constants

def checkOutratesForConvergence(iterations, current_iteration, convergence_range, max_rate_rmsf, debug):
    '''
        checks outgoing rates of all bins for convergence. 
        If all outgoing rates are converged the bin is
        flagged as having converged outrates.
        
        Convergence is achieved if the rmsf of the
        last <convergence_range> iterations is below threshold.
    '''
    iteration_counter = len(iterations) - 1
    # if we can have converged bins (we need <conv_range> iterations. iterations[it_counter] also counts)
    #~ if iterations[iteration_counter].getId() - iterations[0].getId() + 1 >= convergence_range and 
    if len(iterations) >= convergence_range:
        """
           calculate mean rate matrices averaged 
           over the last <reweighting_range>:=rr, rr-1, rr-2 etc
           iterations
        """
        rate_matrices = []
        # calculate mean rate matrices
        for i in range(convergence_range):
            rate_matrices.append(iterations[iteration_counter - (convergence_range-1) + i].RateMatrix())
            #~ print(rate_matrices[-1])
            
        # check bins that exist long enough to have converged rates
        bins_to_check = iterations[iteration_counter - \
        							  (convergence_range-1)].bins[:]
        # check for rate-converged bins        
        for _bin in bins_to_check:
            bin_index = _bin.getId()
            if not iterations[-1].bins[bin_index].isConverged():
                # read rates
                # all_bin_rates is an array of rates to all target_bins
                all_bin_rates = []
                # get rates to all bins that are in the latest matrix.
                for target_bin_index in range(len(rate_matrices[-1][bin_index])):
                    _bin_outrates = []
                    for m in rate_matrices:
                        if len(m[bin_index]) > target_bin_index:
                            _bin_outrates.append(m[bin_index][target_bin_index])
                        else:
                            _bin_outrates.append(0.0)
                    all_bin_rates.append(np.array(_bin_outrates))
               
                # check convergence
                b_converged = True
                neglected_rates = 0
                non_zero_rates = 0    
                
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
                            sys.stderr.write("relative_rmsf({}) > num_boundary but mean{} < num_boundary!. "\
                                            "Interpreted as not converged.\n".format(rmsf, mean))
                    else:
                        relative_rmsf = rmsf/mean
                    
                    if relative_rmsf >= max_rate_rmsf:
                        if mean < 0.5/_bin.getTargetNumberOfSegments():
                            neglected_rates += 1
                        else:
                            b_converged = False
                            break
                    if mean > constants.num_boundary:
                        non_zero_rates += 1
                        
                if b_converged:
                    current_iteration.bins[bin_index].setConverged(True)
                    print ("Bin {bin_id} is converged")
                    if debug:
                        print("{nconv} rates from {nrates} converged. rest neglectable)".format(bin_id = bin_index, nconv=non_zero_rates-neglected_rates, nrates=non_zero_rates))
                					
        #~ print ("rates[0]:")				
        #~ print (rates[0])
    return
