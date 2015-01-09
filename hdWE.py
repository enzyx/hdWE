#!/usr/bin/python3

import argparse

#############################
# parsing the commandline
parser = argparse.ArgumentParser(description=
    'hdWE is a hyperdimensional weighted ensemble simulation implementation.')
parser.add_argument('-d', '--dir', type=str, 
                    dest="workdir", required=True, metavar="DIR",
                    help="The working direcory")
parser.add_argument('-t', '--top', type=str, dest="input_md_topol", 
                    required=True, metavar="FILE",
                    help="System topology file")
parser.add_argument('-p', '--parm', type=str, dest="input_md_param", 
                    required=True, metavar="FILE",
                    help="MD paramter file")
parser.add_argument('-c', '--conf', type=str, dest="input_md_conf", 
                    required=True, metavar="FILE",
                    help="The starting structure file")
parser.add_argument('--trajectories-per-bin', type=int, dest="input_traj_per_bin", 
                    metavar="10", const=10, nargs='?',
                    help="Number of trajectories per bin")
parser.add_argument('--iterations', type=int, dest="input_iterations", 
                    metavar="50", const=50, nargs='?',
                    help="Number of iterations")
parser.add_argument('--threshold', type=float, dest="input_rmsd_threshold", 
                    metavar="0.1", const=0.1, nargs='?',
                    help="Defines the minimal RMSD of a trajectory to all other bins "
                         "after which a new bin is created")
parser.add_argument('--minimal-probability', type=float, dest="input_minimal_probability", 
                    metavar="0.01", const=0.01, nargs='?',
                    help="Minimal probability a trajectory must have to"
                    " allow forking a new bin")
args = parser.parse_args()
#############################

print(args.workdir)
print(args.input_md_topol)
print(args.input_md_param)
print(args.input_md_conf)
