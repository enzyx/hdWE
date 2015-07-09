import os
import shutil
import glob
from lib.iteration import Iteration
from lib.segment import Segment

def prepare(WORKDIR, JOBNAME, starting_structure, OVERWRITE, APPEND, DEBUG):
    """
    Creates the directory structure. Copies the starting configuration
    into the bin_refcoords folder as the first bin.
    """
    sub_dirs = ['{}-run'.format(JOBNAME), '{}-log'.format(JOBNAME), '{}-debug'.format(JOBNAME)]
    # Go to Workdir
    try:
        os.chdir(WORKDIR)
    except OSError:
        raise Exception("Workdir ({}) does not exist.".format(WORKDIR))    
        
    # Create directories if necessary 
    if APPEND == False:
        if not OVERWRITE:
            # Get the backup index 
            backup_index = -1
            for sub_dir in sub_dirs:
                dir_tmp = WORKDIR + sub_dir
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
            
            # Now move old data to backup
            for sub_dir in sub_dirs:               
                dir_tmp = WORKDIR + sub_dir
                try:
                    if backup_index >= 0:
                        os.rename(dir_tmp, "{0}.bak.{1}".format(dir_tmp, backup_index))
                except OSError:
                    print "Could not backup directory {0}".format(dir_tmp)
                    
        # delete old simulation data if overwrite flag is set
        if OVERWRITE:
            for sub_dir in sub_dirs:
                dir_tmp = WORKDIR + sub_dir
                if os.path.exists(dir_tmp):
                    shutil.rmtree(dir_tmp, ignore_errors=True)
                
        # setup new folders and startfile
        for sub_dir in sub_dirs:
            if 'debug' in sub_dir and not DEBUG:
                continue
            os.mkdir(WORKDIR + sub_dir)
        shutil.copyfile(WORKDIR + starting_structure, "{wd}/{jn}-run/00000_00000_00000.rst7".format(wd=WORKDIR, jn=JOBNAME))
        
def parseInitialBoundaries(config):
    initial_boundaries = []
    config_string = config.get('hdWE', 'boundaries')
    for boundary_string in config_string.split(','):
        initial_boundaries.append([])
        for value in boundary_string.split():
            initial_boundaries[-1].append(float(value))
    return initial_boundaries

def parseState(config_string):
    end_state = []
    for coordinate_string in config_string.split('|'):
        end_state.append([])
        for coordinate_ids in coordinate_string.split(','):
            end_state[-1].append(int(coordinate_ids.strip()))        
    
    return end_state

def createInitialIteration(target_number_of_segments, initial_boundaries, md_module):
    """
    Creates the first bin, with the starting structure as refcoords
    and n_segs_per_bin trajectories with probability 1/n_segs_per_bin.
    """
    initial_segment = Segment(probability = 0, 
                              parent_iteration_id = 0, 
                              parent_bin_id = 0, 
                              parent_segment_id = 0,
                              iteration_id = 0, 
                              bin_id = 0, 
                              segment_id = 0)
    initial_coordinates = md_module.calcSegmentCoordinates(initial_segment)
    initial_coordinate_ids = initial_segment.getCoordinateIds(initial_boundaries)
    iteration0 = Iteration(iteration_id = 0, boundaries = initial_boundaries)
    iteration0.generateBin(reference_iteration_id    = iteration0.getId(),
                           reference_bin_id          = 0,
                           reference_segment_id      = 0,
                           target_number_of_segments = target_number_of_segments,
                           coordinate_ids            = initial_coordinate_ids)
    iteration0.bins[0].generateSegment(probability         = 1.0,
                                       parent_iteration_id = 0,
                                       parent_bin_id       = 0, 
                                       parent_segment_id   = 0)
    iteration0.bins[0].segments[0].setCoordinates(initial_coordinates)
    iteration0.bins[0].backupInitialSegments()
    return iteration0
