#!/usr/bin/python2
import argparse
import sys
import ConfigParser
import matplotlib.pyplot as plt         ## plot library
import numpy as np
import constants
from iteration import Iteration
from segment import Segment
from logger import Logger
import convergenceCheck
from hdWE_parameters import HdWEParameters


#### classes ####

class binData(object):
    """
        a container for data of a bin for a single iteration
    """
    def __init__(self, iteration_index, bin_index):
        self.iteration_index = iteration_index
        self.bin_index = bin_index
        self.means = [0]
        self.rmsfs = [0]
        self.rel_rmsfs = [0]
        
    def getId(self):
        return self.bin_index

    def extractConvRates(self, it_datas, target_bin_index, convergence_range):
        """
            @return an array of the previous
            rates from self.bin_index to 
            target_bin_index within convergence range
        """
        conv_rates = []
        for it_data in it_datas:
            if it_data.getId() >= self.iteration_index - convergence_range \
               and it_data.getId() < self.iteration_index \
               and len(it_data.rate_matrix) > self.bin_index\
               and len(it_data.rate_matrix[self.bin_index]) > target_bin_index  :
                    conv_rates.append( it_data.rate_matrix[self.bin_index][target_bin_index] )
        return np.array(conv_rates)

    def calculateMeans(self, iterationData, max_bins, convergence_range, rmsf_threshold):
        """
            calculates the running means, rmsfs 
            and rmsfs/means for the last <convergence_range> iterations
            for every target bin
        """
        self.means = []
        self.rmsfs = []
        self.rel_rmsfs = []
        if self.iteration_index > convergence_range:
            for target_bin in range(max_bins):
                conv_rates = self.extractConvRates(iterationData, target_bin, convergence_range)
                if len(conv_rates) >= convergence_range:
                    # append the mean of the last <conv_range> rates
                    self.means.append(np.mean( conv_rates ))
                    # rmsf (sqrt, mean, square)    
                    self.rmsfs.append(np.sqrt(np.mean(np.square( conv_rates - self.means[-1] ))))
                    #~ print ("\nrates taken into account, the mean and rmsf:")
                    #~ print (conv_rates)
                    #~ print (self.means[-1])
                    #~ print (self.rmsfs[-1])
                    #relative rmsf
                    if self.means[-1] >= constants.num_boundary:
                        self.rel_rmsfs.append(float(self.rmsfs[-1]/self.means[-1]))
                    else:
                        self.rel_rmsfs.append(0.0)

        self.means = np.array(self.means)                    
        self.rmsfs = np.array(self.rmsfs)        
        self.rel_rmsfs = np.array(self.rel_rmsfs)
        
    def getMeans(self):
        return self.means
        
    def getRMSFs(self):
        return self.rmsfs
        
    def getRelRMSFs(self):
        return self.rel_rmsfs    
            
class iterationData(object):
    """
        statistical data for an iteration
    """
    def __init__(self, iteration):
        self.iteration_id = iteration.getId()
        self.rate_matrix = iteration.RateMatrix()
        # create binData list
        self.bin_data = []
        for i in range(iteration.getNumberOfBins()):
            self.bin_data.append(binData(self.iteration_id, i))
        
    def getId(self):
        return self.iteration_id
    
    def getBinData(self, bin_index):
        if len(self.bin_data) > bin_index:
            return self.bin_data[bin_index]
        else:
            return binData(0,0)
        
    def __iter__(self):
        """
        Defines the class as iterable.
        """
        self._iter_index = -1
        return self    

    def __next__(self):
        """
        Returns the next element of the array self.bins
        """
        self._iter_index += 1
        if self._iter_index >= len(self.bin_data):
            raise StopIteration
        else:
            return self.bin_data[self._iter_index]
            
##### functions #####
def findConvegenceEvents(iterations):
    """
        finds the iteration of each bin when it switched to converged.
    """
    # holds the raw iterations of convergence. 0 means no converged rates
    all_closing_iterations = [0] * iterations[-1].getNumberOfBins()
    for iteration in iterations:
        for _bin in iteration.bins:
            if _bin.isConverged() and all_closing_iterations[_bin.getId()] == 0:
                all_closing_iterations[_bin.getId()] = iteration.getId()
    closed_bin        = np.array([i         for i,iteration in enumerate(all_closing_iterations) if iteration != 0])
    closing_iteration = np.array([iteration for i,iteration in enumerate(all_closing_iterations) if iteration != 0]) 
    
    return (closed_bin, closing_iteration)
   

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=
    'Checks when bins were closed.')
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0, metavar='INT',
                    help="First iteration to use.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1, metavar='INT',
                    help="Last iteration to to use.")  
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    required=True, default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing.")
parser.add_argument('-p', '--plot', dest="plot", 
                    default=False, action="store_true",
                    help="Plot the result.")
#TODO: <suspended> created a mutually exclusive group to have one plotting command. probably not the best idea.
#~ functions = parser.add_mutually_exclusive_group(required=True)        
parser.add_argument('-w', '--when', dest="when", 
                    default=False, action='store_true',
                    help="Show the iteration of convergence of bins.")            
parser.add_argument('-r', '--rates', dest="rates",
                    default=False, action='store_true',
                    help="Re-calculates the rates. Requires a config file (-c).")                    
parser.add_argument('-i', '--bin_index', dest="bin_index", 
                   type=int, default=0, metavar="INT",
                   help="Index of bin for outrate-calculation.")
parser.add_argument('-c', '--convergence_range', dest="conv_range", 
                   type=int, default=5, metavar="INT",
                   help="convergence range. How many previous iterations are taken into account for the mean.")                                                    
parser.add_argument('-t', '--threshold', dest="rmsf_threshold", 
                   type=float, default=1.5, metavar="FLOAT",
                   help="threshold of convergence for relative rmsf.")                                
# Initialize
args = parser.parse_args()
                  
################################

# Get the iterations from logger module
logger = Logger(args.logfile)
iterations = logger.loadIterations(first = args.first_iteration, 
                                    last = args.last_iteration)
logger.close()

if args.when:  
    (closed_bin, closing_iteration) = findConvegenceEvents(iterations)
                
    if len(closed_bin) == 0:
        sys.stderr.write('\033[1m' + "No closed bins" + '\033[0m\n')
        sys.exit()
    sys.stdout.write('\033[1m' + 'closed bins:' + '\033[0m\n')
    sys.stdout.write("bin \t closed at iteration\n")     
    for _bin, _iteration in zip(closed_bin, closing_iteration):
        sys.stdout.write("{} \t {}\n".format(_bin, _iteration))   
    
    ## Plot if required
    if args.plot:
        sys.stdout.write('\033[1m' + 'Plotting data...' + '\033[0m\n')
        f, ax = plt.subplots(1, 1)
        ax.grid()
        ax.set_xlim([0, np.amax(closed_bin) + 1])
        ax.set_ylim([0, iterations[-1].getId()])
        ax.set_xlabel('bin')
        ax.set_ylabel('iteration')
        ax.scatter(closed_bin, closing_iteration, alpha=0.5, label = 'bin closed')
        ax.legend()
        plt.show() 

if args.rates:
    sys.stdout.write('\033[1m' + 'Calculating rates...' + '\033[0m\n')
    
    iteration_datasets = []
    n_bins = iterations[-1].getNumberOfBins()
    
    for iteration in iterations:
        iteration_datasets.append(iterationData(iteration))
        
    for it_data in iteration_datasets:
        if len(it_data.bin_data) > args.bin_index:
            it_data.bin_data[args.bin_index].calculateMeans(iteration_datasets, 
                                                        n_bins, 
                                                        args.conv_range, 
                                                        args.rmsf_threshold)

    if args.plot:
        sys.stdout.write('\033[1m' + 'Plotting data...' + '\033[0m\n')
        f, ax = plt.subplots(1,1)
        ax.set_xlim([0, len(iterations)])
        ax.set_xlabel('iteration')
        ax.set_ylabel('bin')
        ax.grid()
        
        # find out maximum rate for point sizing
        max_mean = 0.0
        for iter_data in iteration_datasets:
            means = iter_data.getBinData(args.bin_index).means
            if len(means) > 0:
                bin_max_mean = np.max(means)
                if bin_max_mean >= max_mean:
                        max_mean = bin_max_mean
                

        # plot of rate convergencies
        for iter_data in iteration_datasets:
            means = iter_data.getBinData(args.bin_index).getRelRMSFs()
            rmsfs = iter_data.getBinData(args.bin_index).getRMSFs()
            rel_rmsfs = iter_data.getBinData(args.bin_index).getRelRMSFs()
            y = np.arange(0,len(rel_rmsfs))
            x = [iter_data.getId() for i in range(len(means))]
            #~ print ("target {}: xrange {}-{}".format(target_index, x[0], x[-1]))


            
            # color it green/red/orange below thres/above thres/at 1/nsteps
            c = [""] * len(rel_rmsfs)        
            for i,rr in enumerate(rel_rmsfs):
                if rr <= args.rmsf_threshold:
                    c[i]  = "g"
                elif rr == np.sqrt(args.conv_range - 1.0):
                    c[i] = "orange"
                else:
                    c[i]  = "r"
            ax.scatter(x, y, s=10*means/max_mean, color=c, label = 'relative rmsf')
            #~ ax.get_yaxis().set_ticks([])
            #~ ax.legend(loc="lower left")            
            
    plt.show() 
        
