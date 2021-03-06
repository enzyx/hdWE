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

    def calculateMeans(self, iterationData, max_bins, CONV_RANGE, CONV_THRES):
        """
            calculates the running means, rmsfs 
            and rmsfs/means for the last <CONV_RANGE> iterations
            for every target bin
        """
        self.means = []
        self.rmsfs = []
        self.rel_rmsfs = []
        if self.iteration_index > CONV_RANGE:
            for target_bin in range(max_bins):
                conv_rates = self.extractConvRates(iterationData, target_bin, CONV_RANGE)
                if len(conv_rates) >= CONV_RANGE:
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
    all_closing_iterations = [0] * N_BINS
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
parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    required=True, default="logfile.log", metavar="DIR",
                    help="The logdir to read.")
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    type=int, default=0, metavar='INT',
                    help="First iteration to use.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1, metavar='INT',
                    help="Last iteration to to use.")  
parser.add_argument('-c', '--conf', type=str, dest="configfile", nargs='?',
                    metavar="FILE", default=False,
                    help="An optional hdWE configuration file")
parser.add_argument('-i', '--bin_index', dest="bin_index", 
                   type=int, default=0, metavar="INT",
                   help="Index of bin for outrate-calculation.") 
parser.add_argument('-p', '--plot', dest="plot", 
                    default=False, action="store_true",
                    help="Plot the result.")   
parser.add_argument('-w', '--when', dest="when", 
                    default=False, action='store_true',
                    help="Just show the iteration of convergence of bins.")            
parser.add_argument('-C', '--converged', dest="convergence",
                    default=False, action='store_true',
                    help="Checks convergence of individual rates of bin.")
parser.add_argument('-r', '--rates', dest="rates",
                    default=False, action='store_true',
                    help="Re-calculates the outrates.")                      
parser.add_argument('-t', '--transitions', dest="transitions",
                    default=False, action='store_true',
                    help="Count trasitions for convergence check.")                              
# Initialize
args = parser.parse_args()
               
################################

# Get the iterations and parameters
print ('loading iterations')
logger = Logger(args.logdir)
iterations = logger.loadIterations(begin = args.first_iteration, 
                                   end   = args.last_iteration)
if args.configfile:
    CONFIGFILE = args.configfile
else:
    CONFIGFILE = logger.loadConfigFile(iterations[0].getId())
config = ConfigParser.ConfigParser()
config.read(CONFIGFILE)  

####### Variables ########

CONV_RANGE = int(config.get('hdWE','convergence-range'))
CONV_THRES = float(config.get('hdWE','convergence-threshold'))
N_BINS     = iterations[-1].getNumberOfBins() 
rate_matrices = []
print ("range={}, threshold={}".format(CONV_RANGE, CONV_THRES))

#################
# CLOSED STATUS #
#################
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
    sys.exit()

###################
# CALCULATE RATES #
###################
if args.rates or args.convergence:
    sys.stdout.write('\033[1m' + 'Calculating rates...' + '\033[0m\n')
    
    # fill iteration and bin data
    iteration_datasets = []
    
    for iteration in iterations:
        iteration_datasets.append(iterationData(iteration))
        rate_matrices.append(iteration_datasets[-1].rate_matrix)
        
    for it_data in iteration_datasets:
        if len(it_data.bin_data) > args.bin_index:
            it_data.bin_data[args.bin_index].calculateMeans(iteration_datasets, 
                                                        N_BINS, 
                                                        CONV_RANGE, 
                                                        CONV_THRES)
          
##############
# PLOT RATES #
##############
if args.rates:
    sys.stdout.write("non-zero mean rates of last iteration\n"\
                     "for bin {_bin}:\n".format(_bin=args.bin_index))
    sys.stdout.write("target bin \t mean\n")
    means = iteration_datasets[-1].getBinData(args.bin_index).getMeans()
    for target_index, mean in enumerate(means):
        if mean > constants.num_boundary and target_index != args.bin_index:
            sys.stdout.write("{0} \t {1:.5f}\n".format(target_index, mean))
            
    # get number of plots to draw
    number_of_plots = 0
    for it_data in iteration_datasets:
        means = it_data.getBinData(args.bin_index).getMeans()
        number_of_plots_for_iteration = 0
        for target_bin_index,mean in enumerate(means):
            if target_bin_index == args.bin_index:
                continue
            if mean > constants.num_boundary:
                number_of_plots_for_iteration += 1
        if number_of_plots_for_iteration > number_of_plots:
            number_of_plots = number_of_plots_for_iteration
    print ('creating {} plots'.format(number_of_plots))
            
    

    if args.plot:
        # plot data generation
        # list of tuples (bin_id, [x], [y])
        plot_data = []
        for target_bin_index, target_bin in enumerate(iterations[-1].bins):
            if target_bin_index == args.bin_index:
                continue
            b_plot_this_bin = False
            x = []
            y = []
            for it_data in iteration_datasets:
                x.append(it_data.getId())
                if len(it_data.getBinData(args.bin_index).means) > target_bin_index:
                    mean = it_data.getBinData(args.bin_index).means[target_bin_index]
                    y.append(mean)
                    if mean > constants.num_boundary:
                        b_plot_this_bin = True
                else:
                    y.append(0.0)
            if b_plot_this_bin:
                plot_data.append((target_bin_index, x, y))
                    
         # actual plotting
        f, axarr = plt.subplots(len(plot_data),1, sharex=True, sharey=False)
        for plot_index, data in enumerate(plot_data):
            print ('plotting plot {}'.format(plot_index))
            target_bin_index = data[0]
            x = data[1]
            y = data[2]
            axarr[plot_index].plot(x,y)
            axarr[plot_index].set_ylabel('bin {}'.format(target_bin_index))
            plot_index += 1
            
        plt.show()



###############
# CONVERGENCE #
###############
if args.convergence:
    sys.stdout.write('\033[1m' + 'checking convergence of individual rates...' + '\033[0m\n')     

    # plot of rate convergences
    if args.plot:
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
     
        for iter_data in iteration_datasets:
            bin_data = iter_data.getBinData(args.bin_index)
            # plot all bins except the outrate bin
            means = [mean for i,mean in enumerate(bin_data.getMeans()) if i != bin_data.getId()]
            rmsfs = [rmsf for i,rmsf in enumerate(bin_data.getRMSFs()) if i != bin_data.getId()]
            rel_rmsfs = [rel_rmsf for i,rel_rmsf in enumerate(bin_data.getRelRMSFs()) if i != bin_data.getId()]
            y = [i for i,mean in enumerate(bin_data.getMeans()) if i != bin_data.getId()]
            #y = np.arange(0,len(iter_data.bin_data))
            x = [iter_data.getId() for i in range(len(y))]
            #~ print ("target {}: xrange {}-{}".format(target_index, x[0], x[-1]))
    
    
            
            # color it green/red/orange below thres/above thres/at 1/nsteps
            c = [""] * len(rel_rmsfs)        
            for i,rr in enumerate(rel_rmsfs):
                if rr <= CONV_THRES:
                    c[i]  = "green"
                elif rr == np.sqrt(CONV_RANGE - 1.0):
                    c[i] = "orange"
                else:
                    c[i]  = "red"
            # check for small rates
            #~ for i,mean in enumerate(means):
                #~ if mean > constants.num_boundary:
                    #~ for iteration in iterations:
                        #~ if iteration.getId() == iter_data.getId():
                            #if mean < iteration.bins[i].getProbability()/iteration.bins[i].getNumberOfSegments():
                            #~ if mean < 0.5/(iteration.bins[i].getTargetNumberOfSegments()*iteration.bins[i].getProbability()):
                                #~ c[i] = "blue"
            #size
            size = [0] * len(means)
            for i,mean in enumerate(means):
                if mean > constants.num_boundary:
                    size[i] = 500*mean/max_mean
              
            ax.scatter(x, y, s=size, color=c, marker="s", label = 'relative rmsf')        
        # run our simulation check:
        sys.stdout.write('\033[1m' + 'running convergenceCheck...' + '\033[0m\n')
    
        # compare to convergence result of method from convergenceCheck.py
        bin_converged = [0.0] * len(iterations)
        for iteration in iterations[CONV_RANGE:]:
            it_counter = iteration.getId()
            sys.stdout.write('checking iteration {}\n'.format(it_counter))
            sys.stdout.flush()
            b_converged = convergenceCheck.checkBin(iterations[-1].bins[args.bin_index],
                                                    rate_matrices[it_counter-CONV_RANGE:it_counter+1],
                                                    CONV_THRES,
                                                    debug=False)
            if len(iteration.bins) > args.bin_index \
               and b_converged:
                bin_converged[it_counter] = 1.0
        sys.stdout.write("\n")
    
        # plot if bin was rated as converged by convergenceCheck
        c_conv = [0.0] * len(iterations)
        for i,convergence in enumerate(bin_converged):
            if convergence == 1.0:
                c_conv[i]  = "green"
            else:
                c_conv[i] = "red"
        x = [ iteration.getId() for iteration in iterations ]
        y = [-5] * len(bin_converged)
    
        ax.scatter(x,y, s=50, color=c_conv, marker='o', label = 'convergence')

        # legend
        mark1 = mlines.Line2D(range(1), range(1), color="white", marker='s', markerfacecolor="green")
        mark2 = mlines.Line2D(range(1), range(1), color="white", marker='s', markerfacecolor="red")
        mark3 = mlines.Line2D(range(1), range(1), color="white", marker='s', markerfacecolor="orange")
        mark4 = mlines.Line2D(range(1), range(1), color="white", marker='s', markerfacecolor="blue")
        mark5 = mlines.Line2D(range(1), range(1), color="white", marker='s', markerfacecolor="green", markersize=2)
        mark6 = mlines.Line2D(range(1), range(1), color="white", marker='s', markerfacecolor="green", markersize=4)
        mark7 = mlines.Line2D(range(1), range(1), color="white", marker='s', markerfacecolor="green", markersize=8)
        mark8 = mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="green", markersize=10)
        mark9 = mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="red", markersize=10)
        #~ ax.legend((mark5,mark6,mark7,mark8,mark9,mark1,mark2,mark3,mark4),
                  #~ ('rate < 0.01', 'rate > 0.1', 'rate > 0.1', 
                   #~ 'bin converged', 'bin not converged', 'rate converged','rate not converged', 'just one rate', 
                   #~ 'low mean rate'),
                  #~ numpoints=1, 
                  #~ ncol = 2,
                  #~ loc="upper left")
        ax.legend((mark1,mark2,mark3,mark4, mark8, mark9),
                  ('rate converged','rate not converged', 'just one rate', 
                   'low mean rate', 'bin converged', 'bin not converged'),
                  numpoints=1, 
                  ncol = 2,
                  loc="upper left")   
        
        sys.stdout.write('\033[1m' + 'Plotting data...' + '\033[0m\n')            
        plt.show() 


#########################
### transition counts ###
#########################

if args.transitions:
    print ('analysing transitions')
    transitions_matrix = convergenceCheck.CalculateTransitionsMatrix(iterations)
      
    # PRINT TRANSITIONS MATRIX        
#     for iline,line in enumerate(transitions_matrix):
#         for iele,element in enumerate(line):
#             if iele != iline:
#                 sys.stdout.write("{0:<5d} ".format(element))
#             else:
#                 sys.stdout.write("    0 ")
#         sys.stdout.write("\n")

    # Print all values of the transition matrix
    sys.stdout.write('#from  to    transitions\n')
    for source_bin,line in enumerate(transitions_matrix):
        for target_bin,element in enumerate(line):
            if source_bin != target_bin and element != 0:
                sys.stdout.write("{0:<3d}   {1:<3d}    {2:<4d}\n".format(source_bin, target_bin, element))
            


    # Print transitions of selected bin args.bin_index
#     total_in = 0
#     total_out = 0
#     sys.stdout.write('# total transitions of bin {} \n'.format(args.bin_index))
#     sys.stdout.write('#bin   incoming     outgoing\n')
#     for bin_id in N_BINS:
#         incoming = transitions_matrix[args.bin_index][bin_id]
#         outgoing = transitions_matrix[bin_id][args.bin_index]
#         total_in += incoming
#         total_out += outgoing
#         if incoming != 0 or outgoing != 0:
#             sys.stdout.write('{bin_id:>3}      {inc:>5}      {out:>5}\n'.format(bin_id = bin_id,
#                                                                                 inc    = incoming,
#                                                                                 out    = outgoing))
#     sys.stdout.write('# total numbers:\n')
#     sys.stdout.write('#         {inc:>5}      {out:>5}'.format(inc = total_in,
#                                                                out = total_out))

    # total in and out transitions (from different bins) of all bins
    # column sum:
#     in_transitions  = transitions_matrix.sum(axis=0) - np.diagonal(transitions_matrix)
#     # row sum:
#     out_transitions = transitions_matrix.sum(axis=1) - np.diagonal(transitions_matrix)
# 
#     # print total in and out transitions 
#     sys.stdout.write('# total incmoing and outgoing segments for all bins\n')
#     sys.stdout.write('#bin_id   in      out      in-out\n')
#     for bin_id in range(N_BINS):
#         sys.stdout.write('{bin:<4d}      {inc:<6d}   {out:<6d}   {net:<6d}\n'.format(bin = bin_id, 
#                                                                       inc = in_transitions[bin_id],
#                                                                       out = out_transitions[bin_id],
#                                                                       net = in_transitions[bin_id]-out_transitions[bin_id]))

    # Number of Connected bins
#     connected_bins = []
#     for bin_id in range(N_BINS):
#         transition_partners = 0
#         for target_bin_id, n_transitions in enumerate(transitions_matrix[bin_id]):
#             if n_transitions != 0 and target_bin_id != bin_id:
#                 transition_partners += 1
#         connected_bins.append(transition_partners)
#         
#     # print number of connected bins
#     sys.stdout.write('#number of bins to which each bin hat a segment transfer\n')
#     sys.stdout.write('#bin   number of target bins')
#     for bin_id,n_connected in enumerate(connected_bins):
#         sys.stdout.write('{bin:<3}   {n:<6d}\n'.format(bin = bin_id,
#                                                    n = n_connected))
            

# # Print number of transitions into bins 
#     in_transitions = sorted(in_transitions,key=lambda x: x[1], reverse=True)
#     sys.stdout.write('#total number of incoming segments\n')
#     for bin_data in in_transitions:
#         sys.stdout.write('{bin:<4d}   {transitions:<6d}\n'.format(bin=bin_data[0], transitions=bin_data[1]))
# 
# # Print number of transitions out of bins 
#     out_transitions = sorted(out_transitions,key=lambda x: x[1], reverse=True)
#     sys.stdout.write('#total number of outgoing segments\n')
#     for bin_data in out_transitions:
#         sys.stdout.write('{bin:<4d}   {transitions:<6d}\n'.format(bin=bin_data[0], transitions=bin_data[1]))

                
# TRANSITIONS OUT OF BIN RESOLVED WITH TARGET BIN
#     transition_lists = [[] for i in range(N_BINS)]
# 
#     for it_id,iteration in enumerate(iterations):
#         for bin_id in range(N_BINS):
#             if len(iteration.bins) > bin_id: 
#                 transitions_in_this_iteration = 0
#                 for segment in iteration.bins[bin_id].initial_segments:
#                     if segment.getParentBinId() == args.bin_index:
#                         if segment.getParentIterationId() == it_id - 1:
#                             rate = segment.getProbability() / \
#                                     iterations[it_id-1].\
#                                      bins[segment.getParentBinId()].\
#                                      getProbability()
#                         else:
#                             sys.stderr.write("WARNING: approximate rate for segment {}.\n".format(segment.getNameString()))
#                             rate = 1.0 / iteration.bins[args.bin_index].getTargetNumberOfSegments()
#  
#                         transitions_in_this_iteration += rate 
#                 if len(transition_lists[bin_id]) == 0:
#                     transition_lists[bin_id].append(transitions_in_this_iteration)
#                 else:
#                     transition_lists[bin_id].append(transition_lists[bin_id][-1] + \
#                                                transitions_in_this_iteration)
#             else:
#                 transition_lists[bin_id].append(0)
#              
#     f, ax = plt.subplots(1,1)
#     for bin_id, bin_transitions in enumerate(transition_lists):
#         if bin_id != args.bin_index and bin_transitions[-1] != 0:
#             x_values = [x for x in range(len(bin_transitions))]
#             ax.plot(x_values,bin_transitions, label = "bin {}".format(bin_id))
#     ax.legend(loc=(1.1,0.5))
#     plt.show()

#     sys.stdout.write('#to bin   transitions\n')
#     for target_bin_id, transitions in enumerate(transitions_matrix[args.bin_index][:]):
#         if transitions != 0:
#             sys.stdout.write("{0:>3}      {1:>5}\n".format(target_bin_id, transitions))

# INCOMING TRANSITIONS PER BIN OVER ITERATIONS
#     transition_lists = [[] for i in range(N_BINS)] 
#     for it_id,iteration in enumerate(iterations):
#         for bin_id in range(N_BINS):
#             if len(iteration.bins) > bin_id:
#                 # count transitions into this bin in this iteration
#                 in_transitions = 0
#                 for initial_segment in iteration.bins[bin_id].initial_segments:
#                     if initial_segment.getParentBinId() != bin_id:
#                         in_transitions += 1
# 
#                 # store number of in_transitions
#                 if len(transition_lists[bin_id]) == 0:
#                     transition_lists[bin_id].append(in_transitions)
#                 else:
#                     transition_lists[bin_id].append(transition_lists[bin_id][-1] + \
#                                                in_transitions)
#             else:
#                 transition_lists[bin_id].append(0)
#               
#     f, ax = plt.subplots(1,1)
#     for bin_id, bin_transitions in enumerate(transition_lists):
#         if bin_id != args.bin_index and bin_transitions[-1] != 0:
#             x_values = [x for x in range(len(bin_transitions))]
#             ax.plot(x_values,bin_transitions, label = "bin {}".format(bin_id))
#     ax.legend(loc=(1.1,0.5))
#     plt.show()