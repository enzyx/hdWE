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

        
def parseInitialBoundaries(config):
    initial_boundaries = []
    config_string = config.get('hdWE', 'boundaries')
    for boundary_string in config_string.split(','):
        initial_boundaries.append([])
        for value in boundary_string.split():
            initial_boundaries[-1].append(float(value))
    return initial_boundaries

def parseSampleRegion(config):
    initial_boundaries = []
    config_string = config.get('hdWE', 'sample-region')
    for boundary_string in config_string.split(','):
        initial_boundaries.append([])
        for value in boundary_string.split():
            initial_boundaries[-1].append(float(value))
    return initial_boundaries

def parseStartingStructures(config):
    """
        can a list of filenames or a wildcard with *
    """
    starting_structures = []
    config_string = config.get('hdWE','starting-structures')
    for start_word in config_string.split():
        starting_structures.extend(glob.glob(start_word))           
    return starting_structures

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
                           md_module):
    """
    Creates the first iteration, with equally probable
    starting structures sorted into bins
    Copies the structure into the run folder
    """
    n_segments = len(STARTING_STRUCTURES)
    # setup first iteration
    iteration0 = Iteration(iteration_id = 0, 
                           boundaries = initial_boundaries,
                           sample_region = initial_sample_region)

    # create bins and segments for starting structures
    for starting_structure in STARTING_STRUCTURES:
        coordinates = md_module.calcCoordinatesOfFile(starting_structure)
        coordinate_ids = bin_classifier.getCoordinateIds(coordinates, initial_boundaries)
        sample_region = iteration0.isInSampleRegion(coordinate_ids)
        if not sample_region:
            n_segments -= 1
            sys.stderr.write("\033[0;31m Warning: \033[0m Given structure is not "\
                             "in sampled region. Won't use file {}\n".format(starting_structure))
            continue

        # find or create appropriate bin
        this_bin = None
        for scan_bin in iteration0:
            if scan_bin.getCoordinateIds() == coordinate_ids:
                this_bin = scan_bin
        if this_bin == None:
            bin_id = iteration0.generateBin(reference_iteration_id    = iteration0.getId(),
                                            reference_bin_id          = iteration0.getNumberOfBins(),
                                            reference_segment_id      = 0,
                                            target_number_of_segments = target_number_of_segments,
                                            coordinate_ids            = coordinate_ids,
                                            sample_region             = sample_region)
            this_bin = iteration0.bins[bin_id]

        
        # create segment
        segment_id = this_bin.generateSegment(probability     = 0.0,
                                            parent_iteration_id = 0,
                                            parent_bin_id       = iteration0.getNumberOfBins(), 
                                            parent_segment_id   = 0)
        this_segment = this_bin.segments[segment_id]
        this_segment.setCoordinates(coordinates)
        namestring = this_segment.getNameString()
               
        # copy structure into run dir
        shutil.copyfile(starting_structure, 
                        "{wd}/{jn}-run/{namestring}.{filetype}".format(\
                                                                        wd         = WORKDIR, 
                                                                        jn         = JOBNAME,
                                                                        namestring = namestring,
                                                                        filetype   = starting_structure.split('.')[-1]))           

    # set correct probabilities and copy to initial bins
    for this_bin in iteration0:
        for this_segment in this_bin:
            this_segment.setProbability(1.0/n_segments)
        this_bin.backupInitialSegments()
    
    return iteration0
