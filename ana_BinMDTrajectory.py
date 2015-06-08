#!/usr/bin/env python2
"""
Sorts a given MD trajectory into Ramachandran bins.
"""
from __future__ import print_function
import numpy
import os
from amber_module import MD_module
import argparse  


###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-c', '--conf', dest="input_md_conf", 
                    type=str, default = 'hdWE.conf',
                    help="MD-Software configuration file")
parser.add_argument('-f', '--file', type=str, nargs='+',
                    dest="trajectories", metavar="TRAJFILE", required=True,
                    help="MD trajectory file")
parser.add_argument('-o', '--output', type=str, dest="outfile", 
                    metavar="OUTFILE", default="rama_bins.dat",
                    help="Rama bins output file")
parser.add_argument('-k', '--keep', dest="keep", 
                    default=False,  action='store_true',
                    help="delete the cpptraj files after processing")
parser.add_argument('-g', '--histogram', dest="histogram", 
                    default='', metavar="HIST",
                    help="Write a histogram of frames per bin to file")

# Initialize
# print('\033[1mAnalyzing Rates.\033[0m')      
args = parser.parse_args()
md_module = MD_module(args.input_md_conf, debug=False)

# Write the cpptraj infile
cpptraj_jobname = os.path.basename(args.trajectories[0])

cpptraj_infile_path = "./{cpptraj}.cpptraj_in".format(cpptraj = cpptraj_jobname)
cpptraj_infile      = open(cpptraj_infile_path, 'w')
for traj in args.trajectories: 
    cpptraj_infile.write('trajin {traj}\n'.format(traj=traj))
cpptraj_outfile_path = "{cpptraj}.cpptraj_out".format(cpptraj = cpptraj_jobname)


for residue in md_module.aa_classifier.getRequiredDihedrals():
    for angle in residue['required_angles']:
        cpptraj_infile.write(md_module.getCpptrajDihedralLine(residue['res_id'], angle, cpptraj_outfile_path))

# Write changes to file
cpptraj_infile.close()

# Run cpptraj
cpptraj_execute_string = ' -p {top} -i {inpath} > /dev/null'.format(
                                                    top=md_module.amber_topology_file, 
                                                    inpath=cpptraj_infile_path)
os.system('cpptraj {execute}'.format(execute=cpptraj_execute_string))
        
# Load cpptraj output as numpy array

dihedral_matrix = numpy.loadtxt(cpptraj_outfile_path) 
# Delete the first entry which refers to the frame index
dihedral_matrix = numpy.delete(dihedral_matrix, 0, axis=1)

rama_bins = []
for dihedrals in dihedral_matrix:
    rama_bins.append(md_module.aa_classifier.getBinId(dihedrals))


if not args.keep:
    os.remove(cpptraj_infile_path)
    os.remove(cpptraj_outfile_path)

if args.histogram != '':
    rama_hist = {}
    old_rama_id = rama_bins[0]
    rama_hist[old_rama_id] = {'count': 1, 'transitions': 0}
    for rama_id in rama_bins[1:]:
        try:
            rama_hist[rama_id]['count'] += 1
        except KeyError:
            rama_hist[rama_id] = {'count': 1, 'transitions': 0}
        # Count transitions
        if old_rama_id != rama_id:
            rama_hist[old_rama_id]['transitions'] += 1
        old_rama_id = rama_id
    # Sort
    rama_hist_sorted = sorted(rama_hist.items(), key=lambda x: x[1]['count'], reverse=True)
    rama_hist_file = open(args.histogram, 'w')
    for i in rama_hist_sorted:
        mean_passage_time = 0.0
        if i[1]['transitions'] > 0:
            mean_passage_time = float(i[1]['count'])/float(i[1]['transitions'])
        rama_hist_file.write("{0:s}   {1: 8d} {2: 8d} {3: 15.8f}\n".format(i[0], 
                                                                i[1]['count'], 
                                                                i[1]['transitions'], 
                                                                mean_passage_time))
    rama_hist_file.close()

# Get the bin id
outfile = open(args.outfile, 'w')
for frame, rama_id in enumerate(rama_bins):
    outfile.write("{0: 8d}\t{1:s}\n".format(frame, rama_id))
outfile.close()