#!/usr/bin/python3
import argparse
from segment import Segment
from logger import Logger
import sys
from thread_container import ThreadContainer
import threading

#~ import numpy as np                      ## numeric library
#~ from scipy.optimize import curve_fit    ## fitting library
#~ import matplotlib.pyplot as plt         ## plot library

#### classes ####
class Trajectory():
    """
    a container for complete trajectories
    """
    def __init__(self, segment):
        self.last_segment = Segment(probability         = segment.getProbability(),
                                    parent_iteration_id = segment.getParentIterationId(),
                                    parent_bin_id       = segment.getParentBinId(),
                                    parent_segment_id   = segment.getParentSegmentId(),
                                    iteration_id        = segment.getIterationId(),
                                    bin_id              = segment.getBinId(),
                                    segment_id          = segment.getId())
        self.read_segments = []
        self.read_segments.insert(0,self.last_segment)
    
    def addSegment(self, segment):
        self.read_segments.append(self.last_segment)
        
    def revertSegments(self):
        tmpsegment = self.read_segments[::-1]
        self.segments = tmpsegment
        

##### functions #####

def getIteration(iteration_list, iteration_id):
    """
    searches an iteration in a list with it's id position first.
    if it's not found there all the iteration list is looped over
    @return iteration from iteration list with iteration_id
    """
    if len(iteration_list) > iteration_id \
        and iteration_list[iteration_id].getId() == iteration_id:
        return iteration_list[iteration_id]
    
    for iteration in iteration_list:
        if iteration.getId() == iteration_id:
            return iteration
            
def compileTrajectory(iterations, trajectories, traj_starts, end):
    trajectory = []
    segment = end
    trajectory.append(segment)
    while segment not in traj_starts:
        segment = getIteration(iterations, segment.getParentIterationId())\
                  .bins[segment.getParentBinId()].segments[segment.getParentSegmentId()]
        trajectory.append(segment)
    # add reverted list to trajectories
    trajectories.append(trajectory[::-1])
    return
    
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
    
        

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=
    'Follows segments on their path through bins.')
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

#~ sys.stderr.write('\033[1mReading number of bins per iteration...\033[0m\n')

# Initialize
args = parser.parse_args()

# Get the actual Iteration from logger module
logger = Logger(args.logfile)
iterations = logger.loadIterations(first = args.first_iteration, 
                                    last = args.last_iteration)
logger.close()

thread_container = ThreadContainer()

trajectories = []
print("finished reading data")
        
#~ # find splits and merges
sys.stdout.write("finding starts and ends of trajectories...\n")
traj_starts = []
traj_ends = []
for segment in iterations[0].getSegments():
    traj_starts.append(segment)

prev_seg_list = iterations[0].getSegments()
for iteration in iterations[1:]:    
    """
    go through iterations and see when previous iterations have been 
    used more than once or not at all
    """
    # get lists of previous and current segments
#~ #    #~ for prev_iteration in iterations:
#~ #        #~ if prev_iteration.getId() == iteration.getId()-1:
#~ #            #~ prev_seg_list = prev_iteration.getSegments()
#~ #            #~ break
    seg_list = iteration.getSegments()

    # create list of descendants for each segment
    for prev_segment in prev_seg_list:
        segment_descendant_list=([segment for segment in seg_list if prev_segment.isParent(segment)])
        # check if it's an end
        if len(segment_descendant_list) == 0:
            traj_ends.append(prev_segment)
        # set multiple descendants as start
        traj_starts.extend(segment_descendant_list[1:])
        
    # shifting of segment list:
    prev_seg_list = seg_list

# add final segments as ends
traj_ends.extend(iterations[-1].getSegments())

# print number of starts and ends
print("found starts and ends: \n{start} starts and\n{end} ends".format(\
        start = len(traj_starts),
        end = len(traj_ends)))

# create list of trajectories
sys.stdout.write("put together trajectories...\n")
sys.stdout.write("\n")
trajectories = []
for traj_index,end in enumerate(traj_ends):
    sys.stdout.write("\rcompiling trajectory {}".format(traj_index))
    sys.stdout.flush()
    thread_container.appendJob(threading.Thread(target=compileTrajectory(iterations, trajectories, traj_starts, end)))
    if thread_container.getNumberOfJobs() >= 4:
        thread_container.runJobs()
# Run remaining jobs
thread_container.runJobs()
sys.stdout.write("\n")

# get average length
avg = 0
for trajectory in trajectories:
    avg += len(trajectory)
avg /= len(trajectories)
print("average trajectory length: {avg}".format(avg=avg))

out = open("trajectory_lengths.dat", "w")
for index,traj in enumerate(trajectories):
    out.write("{idx:05d}   {length}\n".format(idx = index, length = len(traj)))
out.close()

    



