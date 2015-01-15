import argparse

#############################
# parsing the commandline
parser = argparse.ArgumentParser(description=
    'hdWE is a hyperdimensional weighted ensemble simulation implementation.')
parser.add_argument('-d', '--dir', type=str, 
                    dest="work_dir", required=True, metavar="DIR",
                    help="The working direcory")
parser.add_argument('-c', '--conf', type=str, dest="input_md_conf", 
                    required=True, metavar="FILE",
                    help="The starting structure file")
parser.add_argument('--segments-per-bin', type=int, dest="segments_per_bin", 
                    metavar="10", default=10, nargs='?',
                    help="Number of trajectories per bin")
parser.add_argument('--iterations', type=int, dest="max_iterations", 
                    metavar="50", default=50, nargs='?',
                    help="Maximum number of iterations")
parser.add_argument('--threshold', type=float, dest="coordinate_threshold", 
                    metavar="0.1", default=0.1, nargs='?',
                    help="Defines the minimal RMSD of a trajectory to all other bins "
                         "after which a new bin is created")
parser.add_argument('--minimal-probability', type=float, dest="input_minimal_probability", 
                    metavar="0.01", default=0.01, nargs='?',
                    help="Minimal probability a trajectory must have to"
                    " allow forking a new bin")
parser.add_argument('--debug', dest="debug", action="store_true",
                    default=False, help="Turn debugging on")
parser_mdgroup = parser.add_mutually_exclusive_group(required=True)
parser_mdgroup.add_argument("--amber", dest="amber", action="store_true",
                    default=False)
parser_mdgroup.add_argument("--gromacs", dest="gromacs", action="store_true",
                    default=False)
args = parser.parse_args()
#############################
