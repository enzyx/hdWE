#!/usr/bin/python3
import argparse
from logger import Logger
import sys
import numpy as np                      ## numeric library
from scipy.optimize import curve_fit    ## fitting library
import matplotlib.pyplot as plt         ## plot library

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=
    'Calculates the PMF along an arbitrary coordinate from the data of a hdWE run.')
parser.add_argument('-b', '--first_it', dest="first_iteration",
                    required=False, type=int, default=0,
                    help="First iteration to use for PMF calculation.")                    
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    required=False, type=int, default=-1,
                    help="Last iteration to to use for PMF calculation.")  
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    required=True, default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")
parser.add_argument('-p', '--plot', dest="plot", 
                    required=False, default=False, action="store_true",
                    help="Plot the fit with mathplotlib")
################################

sys.stderr.write('\033[1mReading number of bins per iteration...\033[0m\n')

# Initialize
args = parser.parse_args()

# Get the actual Iteration from logger module
logger = Logger(args.logfile)
iterations = logger.loadIterations(first = args.first_iteration, 
                                    last = args.last_iteration)
logger.close()

def logistic_curve(t, G, k):
    f0 = 1.0      ## We no the start value is 1 bin
    return G/( 1.0 + np.exp(-k*G*t)*(G/f0 - 1) )

xdata = []
ydata = []

for iteration in iterations:
    xdata.append(iteration.getId())
    ydata.append(iteration.getNumberOfBins())

xdata = np.array(xdata)
ydata = np.array(ydata)

sys.stderr.write('\033[1m' + 'Fitting logistic function...' + '\033[0m\n')

## Fit using startvalue p0 = [A0,B0,C0] and 
## equal uncertainties for data (sigma = None)
p0 = [1.5, 1.0]
popt, pcov = curve_fit(logistic_curve, xdata, ydata, sigma=None, p0=p0, maxfev= 10000)

sys.stderr.write('\033[1m' + 'Fit results:' + '\033[0m\n')
sys.stderr.write('\033[1m' + "G = {0:.0f} ".format(popt[0]) 
                    + '\033[0m' + '(max bins)\n')
sys.stderr.write('\033[1m' + "k = {0:.8f} ".format(popt[1]) 
                    + '\033[0m' + '(proportionality constant)\n')

## Plot if required
if args.plot:
    ## Generate list with x-points for fit plot
    xfit = np.linspace(0, 200, 200)
    sys.stderr.write('\033[1m' + 'Plotting data...' + '\033[0m\n')
    plt.plot(xdata , ydata                    , 'rx', label = 'real')
    plt.plot(xfit, logistic_curve(xfit, *popt), 'b-', label = 'fit')
    plt.legend()
    plt.show()
else:
    sys.stderr.write('\033[1m' + 'Printing data...' + '\033[0m\n')
    for i in range(len(xdata)):
        print(xdata[i], ydata[i])


#for iteration in iterations:
#    print("{0: 5d} {1: 10d}".format(iteration.getId(),
#                                    iteration.getNumberOfBins() ))


