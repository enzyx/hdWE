"""
The config parser sanetizes config strings and input data.
It is not a data structure containing config values.
"""

import glob

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
        reweigthing_max_iteration = int(config.get('hdWE', 'reweigthing-max-iteration'))
    except:
        pass
    return reweigthing_max_iteration   

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