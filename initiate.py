import os
import glob
from iteration import Iteration


def prepare(work_dir, starting_structure, append, debug):
    """
    Creates the directory structure. Copies the starting configuration
    into the bin_refcoords folder as the first bin.
    """
    # Create directories if necessary
    if append == False:
        for sub_dir in ['run', 'log', 'debug']:
            if sub_dir=='debug' and debug==False:
                break
            dir_tmp = work_dir + sub_dir
            if os.path.exists(dir_tmp):
                if glob.glob(dir_tmp + '.bak.*'):
                    backups = glob.glob(dir_tmp + '.bak.*')
                    backup_number = int(sorted(backups)[-1].split(".")[-1]) + 1
                else:
                    backup_number = 1
                os.rename(dir_tmp, dir_tmp + '.bak.' + str(backup_number))
            os.mkdir(dir_tmp)
        os.system('cp ' + work_dir + starting_structure + ' ' + work_dir + 'run/' + '00000_00000_00000.rst7')
        


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
