"""
The config parser sanetizes config strings and input data.
It is not a data structure containing config values.
"""

import glob
import sys

CONFIG_ERROR_CODE = -22

def parseInitialBoundaries(config):
    initial_boundaries = []
    config_string = config.get('hdWE', 'boundaries')
    for boundary_string in config_string.split(','):
        initial_boundaries.append([])
        for value in boundary_string.split():
            initial_boundaries[-1].append(float(value))
    for dimension in initial_boundaries:
        dimension.sort()        
    
    return initial_boundaries

def parseSampleRegion(config):
    sample_region = []
    config_string = config.get('hdWE', 'sample-region')
    for sample_region_string in config_string.split(','):
        sample_region.append([])
        for value in sample_region_string.split():
            sample_region[-1].append(float(value))
    for dimension in sample_region:
        dimension.sort()    
        
    return sample_region

def parseStartBinCoordinateIds(config):
    start_bin_coordinate_ids = []
    config_string = config.get('hdWE', 'start-bin-coordinate-ids')
    for value in config_string.split():
        start_bin_coordinate_ids.append(int(value))
    
    return start_bin_coordinate_ids

def parseStartingStructures(config):
    """
    can a list of filenames or a wildcard with *
    """
    starting_structures = []
    config_string = config.get('hdWE','starting-structures')
    for start_word in config_string.split():
        starting_structures.extend(glob.glob(start_word))           
    return starting_structures

def parseKeepCoordsSegments(config):
    """
    @return number of segments per bin which should be kept
    @default -1: keep all segments
    """
    keep_coords_segments = 0
    try:
        keep_coords_segments = int(config.get('hdWE', 'keep-coords-segments'))
    except:
        pass
    return keep_coords_segments

def parseKeepCoordsFrequency(config):
    """
    @return a integer frequency (in iterations) at which coordinate files are not deleted.
    @default 1: (keep all coordinates)
    """
    keep_coords_frequency = 1
    try:
        keep_coords_frequency = int(config.get('hdWE', 'keep-coords-frequency'))
    except:
        pass
    return keep_coords_frequency

def parseCompressIteration(config):
    """
    @return Compress Iteration flag
    """
    compress_iteration = False
    try:
        compress_iteration = bool(config.get('hdWE', 'compress-iteration').lower() == "true")
    except:
        pass
    return compress_iteration

def parseCompressClosestMask(config):
    """
    When CompressIteration is set, cpptraj can strip away all except 
    those N waters closest to the given mask. See amber manual for additional
    information under (cpptraj, closest)
    @return   Compress closest mask
    @default  Empty mask ""
    """
    compress_closest_mask = ""
    try:
        compress_closest_mask = config.get('hdWE', 'compress-closest-mask').strip()
    except:
        pass
    return compress_closest_mask


def parseReweightingRange(config):
    """
    @return Reweighting range. Default is 0 (reweighting off)
    """
    reweighting_range = 0
    try:
        reweighting_range = float(config.get('hdWE', 'reweighting-range'))
    except:
        pass
    return reweighting_range

def parseReweightingMaxIteration(config):
    """
    @return The iteration until which reweighting is performed. Default is 'inf'
    """    
    reweigthing_max_iteration = 'inf'
    try:
        reweigthing_max_iteration = int(config.get('hdWE', 'reweighting-max-iteration'))
    except:
        pass
    return reweigthing_max_iteration   

def parseNumberOfThreads(config):
    """
    @return  The number of threads which can be used on this machine
    @default 1
    """
    number_of_threads = 1
    try:
        number_of_threads = int(config.get('hdWE','number-of-threads'))
    except:
        pass
    return number_of_threads

def parseMaxIterations(config):
    """
    @return Number of iterations to simulate.
    @default Print error and exit
    """
    try:
        max_iterations = int(config.get('hdWE','max-iterations'))
        return max_iterations
    except:
        print("Error: While processing configuration file entry 'max-iterations'")
        sys.exit(CONFIG_ERROR_CODE)

def parseResamplingMode(config):
    """
    @return {closest, random, weighted, no-merge, split-forward}
    """
    valid_modi = ['closest', 'random', 'weighted', 'no-merge', 'split-forward', 'split-forward-front', 'split-region']
    merge_mode = ""
    try:
        merge_mode = str(config.get('hdWE', 'resampling-mode')).strip()
        if merge_mode not in valid_modi:
            raise BaseException
    except:
        print("Error: Could not find valid resampling mode in config file.")
        print("valid modi are: " + str(valid_modi))
        sys.exit(CONFIG_ERROR_CODE)
    return merge_mode

def parseMergeThreshold(config):
    """
    @return merge threshold
    """
    merge_mode = parseResamplingMode(config)
    merge_threshold = 0
    
    if merge_mode == 'closest':
        try:
            merge_threshold = float(config.get('hdWE', 'merge-threshold'))
        except:
            print("Error: Merge mode 'closest' needs merge-threshold flag in config file.")
    
    return merge_threshold

def parseSegmentPerBin(config):
    """
    @return segments per bin
    """
    try:
        segments_per_bin = int(config.get('hdWE','segments-per-bin'))
        return segments_per_bin
    except:
        print("Error: While processing configuration file entry 'segments-per-bin'.")
        sys.exit(CONFIG_ERROR_CODE)

def parseCalculateVelocities(config):
    """
    @return Boolean whether to calculate velocities for segments
    """
    try:
        calculate_velocity = config.getboolean('hdWE','calculate-velocities')
        return calculate_velocity
    except:
        return False

def parserInitialNumberOfTargetSegments(config):
    """
    @return Target number of segments, a bin should have after its initial creation
    """
    try:
        return config.getint('hdWE','segments-per-bin')
    except:
        print("Error: Need to specify 'segments-per-bin' in config file!")
        sys.exit(CONFIG_ERROR_CODE)
        
def parsePrimaryCoordinate(config):
    """
    @return The coordinate along which the sampling should be increased by splitting
    """
    primary_coordinate = 0
    try:
        primary_coordinate = config.getint('hdWE','primary-coordinate')
    except:
        pass
    return primary_coordinate

def parseSplitRegion(config):
    """
    @return The bin id region in the primary coordinate at which splitting is applied 
    """
    split_region = [0, 9e99]
    try:
        split_region = map(int, (config.get('hdWE','split-region')).strip().split())
    except:
        pass
    return split_region

def parseFrontInterval(config):
    """
    @return Bin range counting back from the first occupied in the front to define front region
    """
    front_interval = 9e99
    try:
        front_interval = config.getint('hdWE','front-interval')
    except:
        pass
    return front_interval
    
def parseSplitForwardNumberOfChildren(config):
    """
    @return number of children a segment is split to, when it moves forward
            in the desired coordinate
    """
    split_forward_number_of_children = 0
    try:
        split_forward_number_of_children = config.getint('hdWE','split-forward-children')
    except:
        pass
    return split_forward_number_of_children

###################
#      AMBER      #
###################

def parseLangevinBinary(config):
    """
    @return  path to amber binary (sander, pmemd, pmemd.cuda)
    @default pmemd (Which is expected to be in $PATH
    """
    langevin_binary = "bbk.py"
    try:
        langevin_binary = config.get('langevin','langevin-binary')
    except:
        pass
    return langevin_binary

def parseLangevinTmpDirPath(config):
    """
    @return  tmp directory path
    @default /tmp 
    """
    tmp_dir = '/tmp'
    try:
        tmp_dir = config.get('langevin', 'tmp-dir')
    except:
        pass
    return tmp_dir

###################
#      AMBER      #
###################

def parseCpptrajBinary(config):
    """
    @return  path to cpptraj binary
    @default cpptraj (Which is expected to be in $PATH
    """
    cpptraj = "cpptraj"
    try:
        cpptraj = config.get('amber','cpptraj-binary')
    except:
        pass
    return cpptraj

def parseAmberBinary(config):
    """
    @return  path to amber binary (sander, pmemd, pmemd.cuda)
    @default pmemd (Which is expected to be in $PATH
    """
    amber_binary = "pmemd"
    try:
        amber_binary = config.get('amber','amber-binary')
    except:
        pass
    return amber_binary

def parseAmberRmsdMask(config):
    """
    @return  rmsd fit mask for closest merge mode
    @default None if mergemode is not closest 
    """
    rmsd_mask = None
    merge_mode = parseResamplingMode(config)
    
    if merge_mode == 'closest':
        try:
            rmsd_mask = str(config.get('amber', 'rmsd-mask'))
        except:
            print("rmsd-mask is required for merge mode 'closest'")
            sys.exit(CONFIG_ERROR_CODE)
    return rmsd_mask

def parseAmberRmsdFitMask(config):
    """
    @return rmsd fit mask if merge mode closest
    """
    rmsd_fit_mask = None
    merge_mode = parseResamplingMode(config)
    
    if merge_mode == 'closest':
        try:
            rmsd_fit_mask = str(config.get('amber', 'rmsd-fit-mask'))
        except:
            print("rmsd-fit-mask is required for merge mode 'closest'")
            sys.exit(CONFIG_ERROR_CODE)
    return rmsd_fit_mask

def parseAmberTmpDirPath(config):
    """
    @return  tmp directory path
    @default /tmp 
    """
    tmp_dir = '/tmp'
    try:
        tmp_dir = config.get('amber', 'tmp-dir')
    except:
        pass
    return tmp_dir

def parseSteadyState(config):
    """
    @return parse steady state flag 
    """
    steady_state = False
    try:
        steady_state = config.getboolean('hdWE','steady-state')#
    except:
        pass
    return steady_state