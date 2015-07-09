from __future__ import print_function
import numpy
import lib.constants as constants
import lib.analysis_operations as analysis_operations

class Reweighting(object):
    '''
    The reweighting class keeps track of the iterations rate matrices
    and alleviates the need of storing all concerned iterations which
    typically uses a lot of memory. 
    '''
    
    def __init__(self, reweighting_range):
        # List of rate matrices. The user needs to call self.storeRateMatrix(iteration)
        # to keep the list of rate matrices synchronized with the iterations
        self.rateMatrices = []
        self.reweighting_range = reweighting_range
    
    def storeRateMatrix(self, iteration):
        '''
        Save a particular rate matrix to the class history
        '''
        self.rateMatrices.append(iteration.RateMatrix())
    
    def calcChiSquare(self, P, K):
        '''
        Calculates the chi squared deviation from the equilibrium condition 
        There are N*(N-1)/2 equilibriums conditions P_i*k_ji = P_j*k_ij
        @param P: Vector of bin probabilities
        @param K: Rate matrix k_ij 
        ''' 
        chi2 = 0.0
        N = len(P)
        for i in range(0,N-1):
            for j in range(i+1,N):
                chi2 += (P[i] * K[i][j] - P[j]*K[j][i])**2
        return chi2/(N)**2

    def meanRateMatrix(self):
        '''
        Calculates the mean values of the bin-to-bin rates over the 
        reweighting range and returns the mean rate Matrix.
        The mean value of outrates of bins that do not exists over the 
        full reweighting range is calculated only over the existing 
        iterations in order to not underestimate them.
        '''        
        # Initialize mean rate matrix
        N = len(self.rateMatrices[-1])
        mean_rate_matrix = numpy.zeros([N,N], float)  
        bin_exists_since = numpy.zeros([N], int)
        end   = len(self.rateMatrices)
        begin = int(end * self.reweighting_range)
        
        # Sum rates
        for index in range(begin, end):
            mean_rate_matrix_tmp = self.rateMatrices[index]
            for i in range(0, len(mean_rate_matrix_tmp)):
                if max(mean_rate_matrix_tmp[i,:]) > constants.num_boundary:
                    bin_exists_since[i] += 1
                    for j in range(0, len(mean_rate_matrix_tmp)):
                        mean_rate_matrix[i,j] += mean_rate_matrix_tmp[i,j]
                    
        # Divide by number of iterations            
        for i in range(0, len(mean_rate_matrix)):
            if bin_exists_since[i] > 0:
                for j in range(0,len(mean_rate_matrix)):
                    mean_rate_matrix[i,j] =  mean_rate_matrix[i,j] / bin_exists_since[i]
                 
        return mean_rate_matrix

    def reweightBinProbabilities(self, iteration):
        '''
        Reweights the bin probabilities according to the steady state equations, using
        the mean rates over the last reweighting_range*len(iteration) iterations. 
        Bins with no outrates to other bins are excluded from the reweighting.
        The total probablity of the reweighted bins is not changed.
        The bin probablities of skipped bins are not changed.
        '''
        
        # 0. Check if all bins are occupied
        for this_bin in iteration:
            if this_bin.getNumberOfSegments() == 0:
                print("\x1b[31m     Found empty bin {}, skipping reweighting\x1b[0m".format(this_bin.getId()))
                return
       
        # 1. get the Rate Matrix, averaged over the last iteration_range iterations
        mean_rate_matrix  = self.meanRateMatrix()
              
        # 2. check for bins with no outrates or only outrate of 1 to itself and add their
        #    indices to delete_bin_index. Add the indices of bins that will be reweighted
        #    to keep_bin_index
        keep_bin_index   = numpy.zeros([0], int)
        delete_bin_index = numpy.zeros([0], int)
        for i in range(0,len(mean_rate_matrix)):
            if max(mean_rate_matrix[i,:]) > 0.0 and max(mean_rate_matrix[i,:]) < 1.0:
                keep_bin_index = numpy.append(keep_bin_index, i)
            else: 
                delete_bin_index = numpy.append(delete_bin_index, i)  
        
        # 3. Reduce mean_rate_matrix by deleting all entries at indices contained in delete_bin_index
        reduced_mean_rate_matrix = numpy.delete(mean_rate_matrix, delete_bin_index, 0)
        reduced_mean_rate_matrix = numpy.delete(reduced_mean_rate_matrix, delete_bin_index, 1)
        
        # 4. Calculate the total Probability of the bins that will be reweighted
        prob_of_reweighted_bins = 0.0
        for i in range(0,len(keep_bin_index)):
            # Nasty fix to handle tuples as probabilities
            if type(iteration.bins[keep_bin_index[i]].getProbability()) == float:
                print('skip')
                continue
            prob_of_reweighted_bins += sum(iteration.bins[keep_bin_index[i]].getProbability())
        
        # 5. Try to solve the steady state equations:
        try:
            bin_probs_from_rates = analysis_operations.BinProbabilitiesFromRates(reduced_mean_rate_matrix)
            # If successfully solved:
            # Rescale the normalized bin probabilities from the steady state equations
            # to the total probability the bins had before:
            for i in range(0,len(keep_bin_index)):
                bin_probs_from_rates[i] = bin_probs_from_rates[i] * prob_of_reweighted_bins
    
            # Skip reweighting if a bin probability would become 0
            # (this could lead to trajectories with 0 probability assigned)
            if numpy.min(bin_probs_from_rates) <= constants.num_boundary:
                print('At least one Bin probability would be zero.')
            else:
                # Assign the new probabilities to the bins. 
                # Keep relative segment probabilities within bins.
                for i in range(0,len(keep_bin_index)):
                    reweight_factor = (1.0 * bin_probs_from_rates[i] / 
                                       sum(iteration.bins[keep_bin_index[i]].getProbability()) )
                    for this_segment in iteration.bins[keep_bin_index[i]]:
                        this_segment.setProbability(this_segment.getProbability() *
                                                    reweight_factor)
    
                # Log: Print the reweighted bin probablities
                #print('     Reweighting: Success. Skipped bins: ' + str(delete_bin_index))
#                 bin_probabilities = iteration.getBinProbabilities()
#                 print('     New chi square (equilibriums condition): {0:g}'.format(
#                     self.calcChiSquare(sum(bin_probabilities), mean_rate_matrix)))
                return
        except:
            # Raise Error if the steady state equations could no successfully be solved
            print('\x1b[31m     Singular rate matrix.\x1b[0m')
        return