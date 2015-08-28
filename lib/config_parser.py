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

def parseMergeMode(config):
    """
    @return {closest, random, weighted, marginonly}
    """
    merge_mode = ""
    try:
        merge_mode = str(config.get('hdWE', 'merge-mode')).strip()
        if merge_mode not in ['closest', 'random', 'weighted', 'none']:
            raise BaseException
    except:
        print("Error: Could not find valid merge mode in config file.")
        sys.exit(CONFIG_ERROR_CODE)
    return merge_mode

def parseMergeThreshold(config):
    """
    @return merge threshold
    """
    merge_mode = parseMergeMode(config)
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
        print("Error: Error: While processing configuration file entry 'segments-per-bin'.")
        sys.exit(CONFIG_ERROR_CODE)
    
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
    merge_mode = parseMergeMode(config)
    
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
    merge_mode = parseMergeMode(config)
    
    if merge_mode == 'closest':
        try:
            rmsd_fit_mask = str(config.get('amber', 'rmsd-fit-mask'))
        except:
            print("rmsd-fit-mask is required for merge mode 'closest'")
            sys.exit(CONFIG_ERROR_CODE)
    return rmsd_fit_mask

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