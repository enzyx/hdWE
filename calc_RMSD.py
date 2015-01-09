import os

def calc_RMSD_AMBER(directory,IT,BIN,SEG,dir_topology):
    """Calculates the RMSD of a configuration corresponding to iteration, bin and segment
    with respect to all bin configurations with the cpptraj of AMBER."""