#!/usr/bin/python2
import argparse
import sys
import ConfigParser
import matplotlib.pyplot as plt         ## plot library
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
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
               and len(it_data.rate_matrix) > self.bin_index:
               
                if len(it_data.rate_matrix[self.bin_index]) > target_bin_index  :   
                    conv_rates.append( it_data.rate_matrix[self.bin_index][target_bin_index] )
                else:
                    # bin wasn't created yet (but could have been, just didn't, so the rate before was 0.0)
                    conv_rates.append( 0.0 )
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
    'Checks when bins were closed or if rates converged. ')
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
                   type=int, nargs='?', default=False, metavar="INT",
                   help="Convergence range. optional, if not set logfile value will be used.")                                                    
parser.add_argument('-t', '--threshold', dest="rmsf_threshold", 
                   type=float, nargs='?', default=False, metavar="FLOAT",
                   help="threshold of convergence for relative rmsf. optional.")                                
# Initialize
args = parser.parse_args()
                  
################################

# Get the iterations and parameters
logger = Logger(args.logfile)
logged_param_string = logger.getHdWEParameterString()
iterations = logger.loadIterations(first = args.first_iteration, 
                                    last = args.last_iteration)
logger.close()
hdWE_parameters = HdWEParameters()
hdWE_parameters.loadJsonParams(logged_param_string)

if args.when:  
    (closed_bin, closing_iteration) = findConvegenceEvents(iterations)
                
    if len(closed_bin) == 0:
        sys.stderr.write('\033[1m' + "No bins converged." + '\033[0m\n')
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
    
    # get parameters
    if args.conv_range:
        conv_range = args.conv_range
    else:
        conv_range = hdWE_parameters.convergence_range
    if args.rmsf_threshold:
        rmsf_threshold = args.rmsf_threshold
    else:
        rmsf_threshold = hdWE_parameters.convergence_threshold

    for iteration in iterations:
        iteration_datasets.append(iterationData(iteration))
        
    for it_data in iteration_datasets:
        if len(it_data.bin_data) > args.bin_index:
            it_data.bin_data[args.bin_index].calculateMeans(iteration_datasets, 
                                                        n_bins, 
                                                        conv_range, 
                                                        rmsf_threshold)
    
    sys.stdout.write("non-zero mean rates of last iteration\n"+
                     "for bin {_bin}:\n".format(_bin=args.bin_index))
    sys.stdout.write("target bin \t mean\n")
    means = iteration_datasets[-1].getBinData(args.bin_index).getMeans()
    for target_index, mean in enumerate(means):
        if mean > constants.num_boundary and target_index != args.bin_index:
            sys.stdout.write("{0} \t {1:.5f}\n".format(target_index, mean))
            

    if args.plot:
        sys.stdout.write('\033[1m' + 'preparing data for plotting...' + '\033[0m\n')
        f, ax = plt.subplots(1,1)
        ax.set_xlim([0, iterations[-1].getId() + 1])
        ax.set_xlabel('iteration')
        ax.set_ylabel('bin')
        #~ ax.grid()
        
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
            bin_data = iter_data.getBinData(args.bin_index)
            # plot all bins except the outrate bin
            means = [mean for i,mean in enumerate(bin_data.getMeans()) if i != bin_data.getId()]
            rmsfs = [rmsf for i,rmsf in enumerate(bin_data.getRMSFs()) if i != bin_data.getId()]
            rel_rmsfs = [rel_rmsf for i,rel_rmsf in enumerate(bin_data.getRelRMSFs()) if i != bin_data.getId()]
            y = [i for i,mean in enumerate(bin_data.getMeans()) if i != bin_data.getId()]
            y = np.arange(0,len(iter_data.bin_data))
            x = [iter_data.getId() for i in range(len(y))]
            #~ print ("target {}: xrange {}-{}".format(target_index, x[0], x[-1]))


            
            # color it green/red/orange below thres/above thres/at 1/nsteps
            c = [""] * len(rel_rmsfs)        
            for i,rr in enumerate(rel_rmsfs):
                if rr <= rmsf_threshold:
                    c[i]  = "green"
                elif rr == np.sqrt(conv_range - 1.0):
                    c[i] = "orange"
                else:
                    c[i]  = "red"
            # check for small rates
            for i,mean in enumerate(means):
                if mean > constants.num_boundary:
                    for iteration in iterations:
                        if iteration.getId() == iter_data.getId():
                            #~ if mean < iteration.bins[i].getProbability()/iteration.bins[i].getNumberOfSegments():
                            if mean < 0.5/iteration.bins[i].getTargetNumberOfSegments():
                                c[i] = "blue"
            #size
            size = [0] * len(means)
            for i,mean in enumerate(means):
                if mean > constants.num_boundary:
                    size[i] = 5
                if mean > 0.01:
                    size[i] = 10
                if mean > 0.1:
                    size[i] = 20
                
            #~ ax.scatter(x, y, s=50*means/max_mean, color=c, marker="o", label = 'relative rmsf')
            ax.scatter(x, y, s=size, color=c, marker="o", label = 'relative rmsf')
            #~ ax.get_yaxis().set_ticks([])
            
        # run our simulation check:
        sys.stdout.write('\033[1m' + 'running convergenceCheck...' + '\033[0m\n')

        bin_converged = [0.0] * len(iterations)
        for it_counter,iteration in enumerate(iterations):
            sys.stdout.write('checking iteration {}\r'.format(it_counter))
            sys.stdout.flush()
            
            #create dummy iteration where convergenceCheck will set its convergence flag
            dummy_iteration = Iteration(iteration.getId())
            parent_iteration = iteration
            # Generate all previous bins for dummy iteration, but set them unconverged
            for parent_bin in parent_iteration:
                dummy_iteration.generateBin(reference_iteration_id=parent_bin.getReferenceIterationId(),
                            reference_bin_id=parent_bin.getReferenceBinId(),
                            reference_segment_id=parent_bin.getReferenceSegmentId(),
                            target_number_of_segments=parent_iteration.getTargetNumberOfSegments(),
                            outrates_converged = False)
            convergenceCheck.checkOutratesForConvergence(iterations[:it_counter+1], dummy_iteration, conv_range, rmsf_threshold)
            if len(dummy_iteration.bins) > args.bin_index \
               and dummy_iteration.bins[args.bin_index].isConverged():
                bin_converged[it_counter] = 1.0

        # color for convergence
        c_conv = [0.0] * len(iterations)
        for i,convergence in enumerate(bin_converged):
            if convergence == 1.0:
                c_conv[i]  = "green"
            else:
                c_conv[i] = "red"
        x = [ iteration.getId() for iteration in iterations ]
        y = [-5] * len(bin_converged)

        ax.scatter(x,y, s=50, color=c_conv, marker='s', label = 'convergence')
        #~ legend_handles = mlines.Line2D(range(1), range(1), color="green", marker='o', label="converged rate")
        #~ line2 = mlines.Line2D(range(1), range(1), color="white", marker='o',markerfacecolor="green")
        #~ ax.legend(handles = legend_handles, loc="upper left")  
        red_patch = mpatches.Patch(color='red', label='The red data')
        #~ plt.legend(handles=[red_patch])
        #~ ax.legend(handles=[red_patch], loc="upper left")
        #~ ax.legend()
        
        #works
        #~ proxy = [plt.Rectangle((0,0),1,1,fc = "r") for i in range(2)]
        #~ plt.legend(proxy, ["range(2-3)", "range(3-4)", "range(4-6)"])
        mark1 = mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="green")
        mark2 = mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="red")
        mark3 = mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="orange")
        mark4 = mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="blue")
        mark5 = mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="green", markersize=2)
        mark6 = mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="green", markersize=4)
        mark7 = mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="green", markersize=8)
        mark8 = mlines.Line2D(range(1), range(1), color="white", marker='s', markerfacecolor="green", markersize=10)
        mark9 = mlines.Line2D(range(1), range(1), color="white", marker='s', markerfacecolor="red", markersize=10)
        ax.legend((mark5,mark6,mark7,mark8,mark9,mark1,mark2,mark3,mark4),
                  ('rate > 0.1', 'rate > 0.01', 'rate < 0.01', 
                   'bin converged', 'bin not converged', 'rate converged','rate not converged', 'just one rate', 
                   'low mean rate'),
                  numpoints=1, 
                  ncol = 2,
                  loc="upper left")
    
        
        sys.stdout.write('\033[1m' + 'Plotting data...' + '\033[0m\n')            
        plt.show() 
        
