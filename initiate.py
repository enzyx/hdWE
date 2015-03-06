import os, sys
import shutil
import glob
from iteration import Iteration


def prepare(work_dir, logfile, starting_structure, overwrite, append, debug):
    """
    Creates the directory structure. Copies the starting configuration
    into the bin_refcoords folder as the first bin.
    """
    sub_dirs = ['run', 'log', 'debug']
    # Create directories if necessary 
    if append == False:
        if not overwrite:
            # Get the backup index 
            backup_index = -1
            for sub_dir in sub_dirs:
                dir_tmp = work_dir + sub_dir
                if os.path.exists(dir_tmp):
                    if backup_index < 0:
                        backup_index = 1
                    backups = glob.glob(dir_tmp + '.bak.*')
                    try:
                        backups = [int(backup.split(".")[-1]) for backup in backups]
                        tmp_index = int(sorted(backups)[-1]) + 1
                        if backup_index < tmp_index:
                            backup_index = tmp_index                  
                    except IndexError:
                        continue
            
            # backup log file
            if os.path.isfile(logfile):
                os.rename(logfile, logfile+".bak." + str(backup_index))
                
            # Now move old data to backup
            for sub_dir in sub_dirs:               
                dir_tmp = work_dir + sub_dir
                try:
                    if backup_index >= 0:
                        os.rename(dir_tmp, "{0}.bak.{1}".format(dir_tmp, backup_index))
                except OSError:
                    print "Could not backup directory {0}".format(dir_tmp)
                    
        # delete old simulation data if overwrite flag is set
        if overwrite:
            for sub_dir in sub_dirs:
                dir_tmp = work_dir + sub_dir
                if os.path.exists(dir_tmp):
                    shutil.rmtree(dir_tmp)
            os.remove(logfile)
                
        # setup new folders and startfile
        for sub_dir in sub_dirs:
            if sub_dir=='debug' and not debug:
                continue
            os.mkdir(work_dir + sub_dir)
        shutil.copyfile(work_dir + starting_structure, work_dir + 'run/' + '00000_00000_00000.rst7')

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
