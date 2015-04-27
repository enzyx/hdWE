#!/usr/bin/env python2
"""
Maps the configurations of a free MD to the bins of a hdWE run.
"""
from __future__ import print_function
import numpy
from logger import Logger
from amber_module import MD_module
import argparse  
import os
import sys

    
###### Parse command line ###### 
parser =argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-c', '--conf', dest="input_md_conf", nargs='?',
                    default=False, 
                    help="Optional MD-Software configuration file")
parser.add_argument('-l', '--log', type=str, dest="logdir", 
                    default="log", metavar="DIR", required=True,
                    help="The log directory")
parser.add_argument('-o', '--output', dest="output_path", 
                    type=str, default='ana_MapToBins.output',
                    help="Output filename")  
parser.add_argument('-t', '--trajectory', dest="trajectory_path", required=True,
                    help="the free MD trajectory.")
parser.add_argument('-e', '--last_it', dest="last_iteration",
                    type=int, default=-1,
                    help="Use bins that exist in this iteration.") 
parser.add_argument('--threshold', dest="threshold",
                    type=int, required=True,
                    help="Threshold used for hdWE bins.") 
                    
# Initialize
print('\033[1mMapping free MD trajectory to hdWE Bins\033[0m')      
args = parser.parse_args()


#get the actual Iteration from logger module
logger = Logger(args.logdir)
iterations = logger.loadIterations(args.last_iteration, args.last_iteration)
last_iteration = iterations[-1]

# load md module
if not args.input_md_conf:
    args.input_md_conf = logger.loadConfigFile(iterations[0].getId())
md_module = MD_module(args.input_md_conf, debug=False)


#write cpptraj file
#load trajectory
cpptraj_infile_path = md_module.workdir+ 'ana_MapToBins.cpptraj_in'
cpptraj_infile      = open(cpptraj_infile_path, 'w')
cpptraj_outfile_path= md_module.workdir + 'ana_MapToBins.cpptraj_out'
cpptraj_infile.write('trajin {trajectory}\n'.format(trajectory=args.trajectory_path))
#calc rmsd to all bin references
for this_bin in last_iteration.bins:
    reference_bin_name_string         = "{jn}-run/{ref_segment}.rst7".format(jn=md_module.jobname,
                                                                         ref_segment=this_bin.getReferenceNameString())
    cpptraj_reference_id_name_string  = '[reference_id_{0:05d}]'.format(this_bin.getId())

    cpptraj_infile.write('reference {0} {1}\n'.format(reference_bin_name_string, cpptraj_reference_id_name_string))
    cpptraj_infile.write('rms in_{0:05d} {1} out {2} ref {3}\n'.format(this_bin.getId(), 
                                                                         md_module.amber_coordinate_mask, cpptraj_outfile_path,
                                                                         cpptraj_reference_id_name_string))
cpptraj_infile.close()

#Run cpptraj
print(' executing cpptraj')
cpptraj_execute_string = ' -p {top} -i {inpath}'.format(
                                                    top=md_module.amber_topology_file, 
                                                    inpath=cpptraj_infile_path)
os.system('cpptraj {execute} > ana_MapToBins.cpptraj_log'.format(execute=cpptraj_execute_string))
        
#Load cpptraj output as numpy array
print(' loading rmsd data')
try:
    data = numpy.loadtxt(cpptraj_outfile_path) 
    #Delete the first entry which refers to the frame index
    data = numpy.delete(data, 0, axis=1)
except:
    #TODO What should happen then?
    print('cpptraj output {0} can not be found or loaded.'.format(cpptraj_outfile_path))

n_frames = len(data[:,0])
#sort configurations into bins
bin_indices = numpy.zeros([len(data[0,:])+1,3], float)
for i in range(0,len(data[:,0])):
    sys.stdout.write(' Sorting frame {frame:05d} /{n_frames}\r'.format(frame = i, n_frames = n_frames))
    sys.stdout.flush()
    frame_handled = False
    for j in range(0,len(data[0,:])):
        #sort into first bin within threshold, bins are chronologically ordered
        if data[i,j] < args.threshold:
            bin_indices[j,0]  += 1
            frame_handled = True
            break
    #if the configuration does not fit into any bin, count into last entry of bin_indices
    if frame_handled == False:
        bin_indices[-1,0] += 1
        
os.remove(cpptraj_infile_path)
os.remove(cpptraj_outfile_path)
#normalize
frame_cnt = 0
for i in range(0,len(bin_indices[:,0])-1):
    frame_cnt += bin_indices[i,0]
for i in range(0,len(bin_indices[:,0])):
    bin_indices[i,1] = 1.0 * bin_indices[i,0] / frame_cnt
for i in range(0,len(bin_indices[:,0])-1):
    bin_indices[i,2] =  last_iteration.bins[i].getProbability()

header = 'free MD frames in bin, normalized to probabiltiy (without bins that didnt fit), hdWE probability per bin'
numpy.savetxt( md_module.workdir + args.output_path, bin_indices )
sys.stdout.write(' Output written to {out} \n'.format(out = md_module.workdir + args.output_path))
sys.stdout.flush()


    