from __future__ import print_function
import sys
import numpy as np
import constants

def checkBin(_bin, rate_matrices, CONV_THRES, debug):
    """
    checks if RMSFs of outgoing rates of 
    last CONV_RANGE iterations of this bin are
    below CONV_THRES. 
    """
    b_converged = True
    bin_index = _bin.getId()
    non_zero_rates = 0
    neglected_rates = 0
    
    # for all bins existent at last given iteration
    for target_bin_index in range(len(rate_matrices[-1][bin_index])):
        # read rates
        outrates = []
        for matrix in rate_matrices:
            if len(matrix) > bin_index and len(matrix[bin_index]) > target_bin_index:
                outrates.append(matrix[bin_index][target_bin_index])
            else:
                outrates.append(0.0)
                
        # do statistics on rates        
        mean = np.mean(outrates)
        rmsf = np.sqrt(np.mean(np.square(outrates - mean)))
        if mean < constants.num_boundary:
            if rmsf < constants.num_boundary:
                relative_rmsf = 0.0
            else:
                # in the exotic case that the mean is basically zero 
                # but the rmsf isn't we don't consider it converged
                relative_rmsf = CONV_THRES + 1
                sys.stderr.write("relative_rmsf({}) > num_boundary but mean{} < num_boundary!. "\
                                "Interpreted as not converged.\n".format(rmsf, mean))
        else:
            relative_rmsf = rmsf/mean

        #print ("{}->{} relative rmsf: {}".format(bin_index, target_bin_index, relative_rmsf))
        if relative_rmsf >= CONV_THRES:
            if mean < 0.5/_bin.getTargetNumberOfSegments():
                neglected_rates += 1
            else:
                b_converged = False
                break
        if mean > constants.num_boundary:
            non_zero_rates += 1
    if debug:
        if neglected_rates >= 1:
            print("Bin {bin_id}: {nconv} rates from {nrates} converged. "\
              "rest neglectable)".format(bin_id = bin_index,
                                         nconv  = non_zero_rates-neglected_rates, 
                                         nrates = non_zero_rates))
        else:
            sys.stdout.write("Bin {}: no neglected rates.\n".format(bin_index))
    return b_converged

def checkOutratesForConvergence(iterations, CONV_RANGE, CONV_THRES, debug):
    '''
        checks outgoing rates of all bins for convergence. 
        If all outgoing rates are converged the bin is
        flagged as having converged outrates.
        
        Convergence is achieved if the rmsf of the
        last <CONV_RANGE> iterations is below threshold.
    '''
    iteration_counter = len(iterations) - 1
    converged_bins = []
    # we need <conv_range> iterations.
    if len(iterations) < CONV_RANGE:
        return
    
    """
       calculate mean rate matrices averaged 
       over the last <reweighting_range>:=rr, rr-1, rr-2 etc
       iterations
    """
    # calculate the rate matrices
    rate_matrices = []
    for i in range(CONV_RANGE):
        rate_matrices.append(iterations[iteration_counter - (CONV_RANGE-1) + i].RateMatrix())
        #~ print(rate_matrices[-1])
        
    # check bins that exist long enough and are not closed
    bins_to_check = []
    for _bin in iterations[iteration_counter - (CONV_RANGE-1)].bins:
        if not iterations[-1].bins[_bin.getId()].isConverged():
            bins_to_check.append(_bin)
            
    # check for rate-converged bins        
    for _bin in bins_to_check:
        bin_index = _bin.getId()
        b_converged = checkBin(_bin, rate_matrices, CONV_THRES, debug)
        if b_converged:
            iterations[-1].bins[bin_index].setConverged(True)
            converged_bins.append(bin_index)
            
    if len(converged_bins) > 0:
        print ("     Newly converged bins: {}".format(converged_bins))        					
    #~ print ("rates[0]:")				
    #~ print (rates[0])
    return
