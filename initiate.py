import os
from iteration import Iteration
from bin import Bin
from segment import Segment

def prepare(work_dir, starting_structure, override, debug):
    """
    Creates the directory structure. Copies the starting configuration
    into the bin_refcoords folder as the first bin.
    """
    # Create or override directories
    # /run    
    try:
        os.mkdir(work_dir + '/run')
        #copy starting_structure to /run directory
    except:
        # if already exists, move run directory to a backup directory 
        if override == True:
            os.rmdir(work_dir + '/run')
    # /debug
    # /log

def create_initial_iteration(target_number_of_segments):
    """
    Creates the first bin, with the starting structure as refcoords
    and n_segs_per_bin trajectories with probability 1/n_segs_per_bin.
    """
    iteration0 = Iteration(iteration_id=0)
    iteration0.generateBin(reference_iteration_id=iteration0.getId(),
                           reference_bin_id=0,
                           reference_segment_id=0,
                           target_number_of_segments=target_number_of_segments)
    iteration0.bins[0].generateSegment(probability=1.0,
                                       parent_bin_id=0, 
                                       parent_segment_id=0)
    
    return iteration0
