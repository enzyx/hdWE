import os, sys
import shutil
import glob
from lib.iteration import Iteration
import lib.bin_classifier as bin_classifier


def prepare(WORKDIR, JOBNAME, OVERWRITE, APPEND, DEBUG):
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
                
        # setup new folders
        for sub_dir in sub_dirs:
            if 'debug' in sub_dir and not DEBUG:
                continue
            os.mkdir(WORKDIR + sub_dir)


def parseState(config_string):
    end_state = []
    for coordinate_string in config_string.split('|'):
        end_state.append([])
        for coordinate_ids in coordinate_string.split(','):
            end_state[-1].append(int(coordinate_ids.strip()))        
    
    return end_state

def createInitialIteration(STARTING_STRUCTURES, 
                           WORKDIR,
                           JOBNAME,
                           target_number_of_segments, 
                           initial_boundaries, 
                           initial_sample_region, 
                           md_module,
                           start_bin_coordinate_ids):
    """
    Creates the first iteration, with equally probable
    starting structures sorted into bins
    Copies the structure into the run folder
    """

    # randomize the order of entries in starting_structures in order to 
    # make the random choice of first starting structures later straight-forward
    from random import shuffle
    shuffle(STARTING_STRUCTURES)

    # setup first iteration
    iteration0 = Iteration(iteration_id = 0, 
                           boundaries = initial_boundaries,
                           sample_region = initial_sample_region,
                           n_starting_structures = 0)

    # create bins and segments for starting structures
    for starting_structure in STARTING_STRUCTURES:
        coordinates    = md_module.calcCoordinatesOfFile(starting_structure)
        coordinate_ids = bin_classifier.getCoordinateIds(coordinates, initial_boundaries)
        sample_region  = iteration0.isInSampleRegion(coordinate_ids)
        # Check whether starting structure is in active start bin
        if not sample_region:
            sys.stderr.write("\033[0;31m Warning: \033[0m Given structure is not "\
                             "in sampled region. Won't use file {}\n".format(starting_structure))
            continue
        if not all(coordinate_ids == start_bin_coordinate_ids):
            sys.stderr.write("\033[0;31m Warning: \033[0m Given structure is not "\
                             "in starting bin. Won't use file {}\n".format(starting_structure))
            continue
        
        # Create start bin from the first structure in starting structures.
        # If the loop arrives here, it is in the starting bin.
        if len(iteration0.bins) == 0:
            bin_id = iteration0.generateBin(target_number_of_segments = target_number_of_segments,
                                            coordinate_ids            = coordinate_ids,
                                            sample_region             = sample_region)
            starting_bin = iteration0.bins[bin_id]
            
        # Create segment from starting structure
        segment_id = starting_bin.generateSegment(probability         = 0.0,
                                                  parent_iteration_id = 0,
                                                  parent_bin_id       = starting_bin.getId(), 
                                                  parent_segment_id   = starting_bin.getNumberOfSegments() )        
        starting_segment = starting_bin.segments[segment_id]
        starting_segment.setCoordinates(coordinates)
        
        # copy structure into run dir
        shutil.copyfile(starting_structure, 
                        "{wd}/{jn}-run/{namestring}.{filetype}".format(\
                                                                        wd         = WORKDIR, 
                                                                        jn         = JOBNAME,
                                                                        namestring = starting_segment.getNameString(),
                                                                        filetype   = starting_structure.split('.')[-1]))           

    iteration0.number_starting_structures = starting_bin.getNumberOfSegments()
    # All appropriate starting structures have been copied into the rundir
    # in iteration 0, bin 0 and segment id (in ascending order)
    # Corresponding segments have been generated in the start bin.
    # delete segments from bin in order to end up with target number of segments
    # as the starting_structure list was randomized, simply all segments with 
    # segment_id >= target number of segments can be deleted
    del starting_bin.segments[target_number_of_segments:]
      
    # set correct probabilities and copy to initial bins
    n_segments = starting_bin.getNumberOfSegments()
    for this_segment in starting_bin.segments:
        this_segment.setProbability(1.0/n_segments)
    starting_bin.backupInitialSegments()
    
    return iteration0
