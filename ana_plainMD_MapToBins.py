#!/usr/bin/python2.7
import argparse

# Parse command line
parser = argparse.ArgumentParser(description=
    'Calculates the PMF based on a plain MD trajectory.')
parser.add_argument('-c', '--conf', dest="input_md_conf", 
                    type=str,
                    help="MD-Software configuration file")
parser.add_argument('-l', '--log', type=str, dest="logfile", 
                    default="logfile.log", metavar="FILE",
                    help="The logfile for reading and writing")
parser.add_argument('-b', '--begin_frame', dest="begin_frame",
                    required=False, type=int, default=0,
                    help="First frame to use for PMF calculation.")                    
parser.add_argument('-e', '--end_frame', dest="end_frame",
                    required=False, type=int, default=0,
                    help="Last frame to to use for PMF calculation.")  
parser.add_argument('-o', '--output', dest="output_path", 
                    type=str, default='ana_calculatePMF.output',
                    help="Output filename")  
parser.add_argument('-N', '--number_of_bins', dest="number_of_bins",
                    type=int, default=100, 
                    help="Number of bins used to calculate the probability histogram.")  
parser.add_argument('-i', '--cpptraj_lines_file', dest="cpptraj_lines_file_path", 
                    type=str, 
                    help="File containig cpptraj syntax that defines the reaction coordinate.")
                    
# Initialize
args = parser.parse_args()
