#!/usr/bin/python2
"""
hdWE is a hyperdimensional weighted ensemble simulation implementation.
"""
import sys
import os
import glob
import numpy
import initiate
import analysis_operations
from iteration import Iteration
from logger import Logger
import threading
from thread_container import ThreadContainer
from hdWE_parameters import HdWEParameters

# Compatibility mode for python2.6
has_argparse = False
try:
    import argparse  
    has_argparse = True  
except ImportError:
    import optparse  #Python 2.6

###### Parse command line ###### 
if has_argparse:
    parser =argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
else:
    parser = optparse.OptionParser()
    parser.add_argument = parser.add_option

parser.add_argument('-d', '--dir', type=str, 
                    dest="work_dir", required=False, metavar="DIR",
                    help="The working direcory")
parser.add_argument('-c', '--conf', type=str, dest="md_conf", 
                    required=False, metavar="FILE",
                    help="The starting structure file")
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")          
parser.add_argument("--append", dest="append", action='store_true', 
                    default=False,
                    help="continue previous iterations from logfile with new parameters.") 
parser.add_argument("--resume", dest="resume", action='store_true', 
                    default=False,
                    help="resume previous run with parameters from logfile.")                                                                      
parser.add_argument('--segments-per-bin', type=int, dest="segments_per_bin", 
                    metavar="50", default=50, nargs='?',
                    help="Number of trajectories per bin")
parser.add_argument('--iterations', type=int, dest="max_iterations", 
                    metavar="50", default=False, const=50, nargs='?',
                    help="Maximum number of iterations")
parser.add_argument('--threshold', type=float, dest="coordinate_threshold", 
                    metavar="0.1", default=0.1, nargs='?',
                    help="Defines the minimal RMSD of a trajectory to all other bins "
                         "after which a new bin is created")
parser.add_argument('--minimal-probability', type=float, dest="minimal_probability", 
                    metavar="0.0", default=0.0, nargs='?',
                    help="Minimal probability a trajectory must have to"
                    " allow forking a new bin")
parser.add_argument('--max-bins', type=int, dest="max_bins", 
                    metavar="N", nargs='?', default=-1,
                    help="Maximum number of bins")                    
parser.add_argument('--debug', dest="debug", action="store_true",
                    default=False, help="Turn debugging on")
parser.add_argument('-nt', '--number-of-threads', type=int, dest="number_of_threads", 
                    metavar="1", default=1, nargs='?',
                    help="Number of threads for script parallelization")
parser_mdgroup = parser.add_mutually_exclusive_group(required=False)
parser_mdgroup.add_argument("--amber", dest="amber", action="store_true",
                    default=False)
parser_mdgroup.add_argument("--gromacs", dest="gromacs", action="store_true",
                    default=False)
args = parser.parse_args()
# guarantee a working work_dir variable
if args.work_dir and args.work_dir[-1] != "/":
    args.work_dir +="/"
# define md_package string (amber/gromacs)
if args.amber:
    md_package = "amber"
elif args.gromacs:
    md_package = "gromacs"
else:
    raise ValueError("No MD package selected!")
# backup log file
if not (args.append or args.resume):
    if os.path.isfile(args.logfile):
        if os.path.isfile(args.logfile + '.bak*'):
            backups = glob.glob(args.logfile + '.bak.*')
            backup_number = int(sorted(backups)[-1].split(".")[-1]) + 1
        else:
            backup_number = 1
        os.rename(args.logfile, args.logfile+".bak." + str(backup_number))
        
#############################
# Functions
#############################
def run_parallel_jobs(job_list):
    """
    Run a list of jobs in parallel
    @param job_list of jobs to run
    @return a new empty job_list
    """
    for job in parallel_jobs:
       job.start()
    # Wait until threads are finished
    for job in parallel_jobs:
        job.join()
    # Reset the job list to fill it with next bunch of work
    return []

#############################
# Main
#############################

# The global list of arrays
iterations = [] 

# Initialize the paramters container
hdWE_parameters = HdWEParameters()
if not (args.append or args.resume):
    hdWE_parameters.loadParams(work_dir             = args.work_dir,
                               md_conf              = args.md_conf,
                               md_package           = md_package,
                               max_iterations       = args.max_iterations,
                               segments_per_bin     = args.segments_per_bin,
                               minimal_probability  = args.minimal_probability,
                               coordinate_threshold = args.coordinate_threshold,
                               max_bins             = args.max_bins,
                               logfile              = args.logfile,
                               debug                = args.debug)

# Initialize the logger
logger = Logger(filename=args.logfile, debug=args.debug)
# append or resume
if args.append:
    logger.logParameters(hdWE_parameters)
    iterations = logger.loadIterations() 
elif args.resume:
    hdWE_parameters.loadJsonParams(logger.getHdWEParameterString())
    if args.max_iterations:
        hdWE_parameters.max_iterations = args.max_iterations
        logger.logParameters(hdWE_parameters)
    iterations = logger.loadIterations()   
else:
    logger.logParameters(hdWE_parameters)
    

# Setup the work_dir and initiate iterations
initiate.prepare(hdWE_parameters.work_dir, starting_structure="", override="", debug=hdWE_parameters.debug)
if len(iterations)==0:
    iterations.append(initiate.create_initial_iteration(hdWE_parameters.segments_per_bin))
    logger.logIteration(iterations[0])

# Check MD suite
if(hdWE_parameters.md_package.lower() == "amber"):
    from amber_module import MD_module
    md_module = MD_module(hdWE_parameters.work_dir, hdWE_parameters.md_conf, debug=args.debug)
if(hdWE_parameters.md_package.lower() == "gromacs"):
    print("Sorry, support for gromacs is not implemented yet.")
    sys.exit(-1)


# Loop
for iteration_counter in range(len(iterations), hdWE_parameters.max_iterations + 1):
    iteration = Iteration(iteration_counter)
    parent_iteration = iterations[iteration_counter - 1]
    # Generate all previous bins for new iteration
    for parent_bin in parent_iteration:
        iteration.generateBin(reference_iteration_id=parent_bin.getReferenceIterationId(),
                              reference_bin_id=parent_bin.getReferenceBinId(),
                              reference_segment_id=parent_bin.getReferenceSegmentId(),
                              target_number_of_segments=hdWE_parameters.segments_per_bin)
    
    coordinates = numpy.array([]) # numpy array?
    # Sort segments in bins. Generate new bin if required
    for parent_bin in parent_iteration:
        for segment in parent_bin:
            coordinates = md_module.CalculateCoordinate(segment, iteration.bins)
            min_coordinate = numpy.min(coordinates)
            # Sort Segment into appropriate bin
            if (min_coordinate <= hdWE_parameters.coordinate_threshold or segment.getProbability() < hdWE_parameters.minimal_probability):
                bin_id = numpy.argmin(coordinates)
                iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                  parent_bin_id=segment.getBinId(),
                                                  parent_segment_id=segment.getId())
            # If necessary create new bin
            else:
                bin_id = iteration.generateBin(reference_iteration_id=segment.getIterationId(),
                                      reference_bin_id=segment.getBinId(),
                                      reference_segment_id=segment.getId(),
                                      target_number_of_segments=hdWE_parameters.segments_per_bin)
                iteration.bins[bin_id].generateSegment(probability=segment.getProbability(),
                                                  parent_bin_id=segment.getBinId(),
                                                  parent_segment_id=segment.getId())
    # Split and merge (Manu)
    # Parallel
    if args.number_of_threads > 1:
        thread_container = ThreadContainer()
        for this_bin in iteration:
            thread_container.appendJob(threading.Thread(target=this_bin.resampleSegments))
            if thread_container.getNumberOfJobs() >= args.number_of_threads:
                thread_container.runJobs()
        # Run remaining jobs
        thread_container.runJobs()
    # Serial
    else:
        for this_bin in iteration:
            this_bin.resampleSegments()

    if hdWE_parameters.debug: 
        print("The overall probabiliy is at {0:05f}".format(iteration.getProbability()))

    # Run MD
    md_module.RunMDs(iteration)
    
    iterations.append(iteration)
    
    #test: print rate matrix, bin probabilities from rates and actual bin probabilities
    if iteration_counter > 15:
        mean_rate_matrix= analysis_operations.meanRateMatrix(iterations, iteration_counter -10, iteration_counter)
        print('')
        try:
            print('Bin probabilities from Rate Matrix:')
            print(analysis_operations.BinProbabilitiesFromRates(mean_rate_matrix)) 
        except:

            print('Singular rate matrix')
        print('Actual Bin probabilities:')
        print(analysis_operations.meanBinProbabilities(iterations, iteration_counter -10, iteration_counter))         

        
    # log iteration (Rainer)
    logger.logIteration(iteration)
    
logger.close()

#count total n of segments
n_segments = 0
for iteration_loop in iterations:
     n_segments += iteration_loop.getNumberOfSegments()
    
print('hdWE completed. Total number of propagated segments: ' + str(n_segments) + '            ')
