from __future__ import print_function
import sys
import numpy as np
import constants

def checkOutratesForConvergence(iterations, current_iteration, convergence_range, max_rate_rmsf):
    '''
        checks outgoing rates of all bins for convergence. 
        If all outgoing rates are converged the bin is
        flagged as having converged outrates.
        
        Convergence is achieved if the rmsf of the
        last <convergence_range> iterations is below threshold.
    '''
    iteration_counter = len(iterations) - 1
    # if we can have converged bins
    if iterations[iteration_counter].getId() - iterations[0].getId() >= convergence_range:
        """
           calculate mean rate matrices averaged 
           over the last <reweighting_range>:=rr, rr-1, rr-2 etc
           iterations
        """
        rate_matrices = []
        # calculate mean rate matrices
        for i in range(convergence_range):
            rate_matrices.append(iterations[iteration_counter - (convergence_range - 1) + i].RateMatrix())
            #~ print(rate_matrices[-1])
            
        # bins that exist long enough to have converged rates
        bins_to_check = iterations[iteration_counter - \
        							  convergence_range].bins[:]
        #check for rate-converged bins        
        # for every bin that exists long enough to have converged rates
        for _bin in bins_to_check:
            bin_index = _bin.getId()
            if not iterations[-1].bins[bin_index].isConverged():
                # check if all rates have existed for <convergence_range>
                b_all_rates_exist = True
                # rates is an array of rates to all target_bins
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
               
                # check convergence
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
                    
                    #~ sys.stdout.write("{0:02d}->{1:02d} mean: {2:.3f} rel_rmsf: {3:.3f}\n".format(
                        #~ bin_index, 
                        #~ target_bin_index,
                        #~ mean,
                        #~ relative_rmsf))
                    if relative_rmsf >= max_rate_rmsf and mean > 0.5/_bin.getTargetNumberOfSegments() :
                        b_converged = False
                        break
                        
                if b_converged:
                    current_iteration.bins[bin_index].setConverged(True)
                    print ("Bin {} is converged\n".format(bin_index))
                					
        #~ print ("rates[0]:")				
        #~ print (rates[0])
    return
