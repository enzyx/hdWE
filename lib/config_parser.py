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