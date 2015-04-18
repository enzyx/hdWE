#!/usr/bin/python2
from logger import Logger
import argparse

# Parse command line
parser = argparse.ArgumentParser(description=
    'Print the filenames of all bin references')
parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    metavar="DIRECTORY", required=True,
                    help="The log directory with simulation data")
parser.add_argument('-n', dest='newline', action='store_true',
                    default=False, help="print filenames without newline")

# Initialize
args = parser.parse_args()

logger = Logger(args.logdir)

iteration = logger.loadIteration(-1)
for this_bin in iteration.bins:
    if args.newline:
        print this_bin.getReferenceNameString() + ".rst7",
    else:
        print this_bin.getReferenceNameString() + ".rst7"

