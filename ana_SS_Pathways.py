#!/usr/bin/python2
import argparse
import sys
import numpy
from lib.logger import Logger
import os
import shutil

class Pathway():
    def __init__(self, segment):
        self.frames = [[segment.getIterationId(), segment.getBinId(), segment.getId()]]
        self.frame_positions = []
        
    def addPreviousFrame(self, segment):
        self.frames.insert(0, [segment.getParentIterationId(), segment.getParentBinId(), segment.getParentSegmentId()])  
    
    def addFramePosition(self, position):
        self.frame_positions.insert(0, position)
        
def calcIterationOfDivergence(pathway_I, pathway_J):
    l = len(pathway_I.frames)
    if pathway_I == pathway_J:
        return -1
    for iteration in range(0, l):
        segment_tmp_I =  pathway_I.frames[iteration]
        segment_tmp_J =  pathway_J.frames[iteration]
        if not (segment_tmp_I[1] == segment_tmp_J[1] and segment_tmp_I[2] == segment_tmp_J[2]):
            return iteration
    return -2

###### Parse command line ###### 
parser = argparse.ArgumentParser(description=
    'Extract continuous from hdWE simulation.')
parser.add_argument('-l', '--log', type=str, dest="logdir",
                    required=True, default="hdWE-log", metavar="DIR",
                    help="The logdir to load.")
parser.add_argument('-r', '--run', type=str, dest="rundir",
                    required=False, default="hdWE-run", metavar="DIR",
                    help="The rundir to load.")
parser.add_argument('-e', '--last_it', dest="last_iteration_id",
                    type=int, default=-1, metavar='INT',
                    help="Last iteration to to read.")
parser.add_argument('-t', '--target-bin', dest="target_bin", metavar="INT",
                    required=True, type=int, 
                    help="Extract pathways arriving in this bin.")  

# Initialize
args              = parser.parse_args()
rundir            = args.rundir.strip('/')
logdir            = args.logdir.strip('/')
target_bin        = args.target_bin
logger = Logger(args.logdir)
last_iteration_id = args.last_iteration_id
if last_iteration_id < 0:
    last_iteration_id = logger.getLastIterationId()
last_iteration = logger.loadIteration(last_iteration_id)
pathways = []

sys.stdout.write('\033[1mana_SS_Pathways\033[0m\n')
sys.stdout.write('  Analyzing\n')
sys.stdout.write('   - Target bin: {} (Coordinate Ids: {})\n'.format(target_bin, 
                                                last_iteration.bins[target_bin].getCoordinateIds()))
sys.stdout.flush()

# Initialize pathways from segments in target bin
for this_segment in last_iteration.bins[target_bin].segments:
    pathways.append(Pathway(this_segment))
sys.stdout.write('   - Found {} pathways ending in target bin in iteration {}\n'.format(len(pathways), last_iteration_id))
sys.stdout.flush()

# Reconstruct pathways from log 
for iteration_id in range(last_iteration_id, -1, -1):
    sys.stdout.write('  Backtracing pathways: Iteration {:05d}\r'.format(iteration_id))
    sys.stdout.flush()
    current_iteration     = logger.loadIteration(iteration_id)
    current_segment_list  = open('{rundir}/{iteration}.segment_list'.format(
                                                     rundir    = rundir, 
                                                     iteration = current_iteration.getNameString())
                                                     ,'r').readlines()
    for pathway in pathways:
        segment_tmp          = current_iteration.bins[pathway.frames[0][1]].segments[pathway.frames[0][2]]
        segment_tmp_name_str = segment_tmp.getNameString()
        index_tmp            = current_segment_list.index(segment_tmp_name_str+'\n')   

        # parent segment of iteration 1 corresponds to iteration 0: already arrived at origin of pathway
        if segment_tmp.getIterationId() > 0:
            pathway.addPreviousFrame(segment_tmp)        
        pathway.addFramePosition([iteration_id, index_tmp])



# Write cpptraj_in files into folder
dirname = 'ana_Pathways_{rundir}_{target_bin:05d}/'.format(rundir = rundir, target_bin = target_bin)
if os.path.exists(dirname):
    shutil.rmtree(dirname)
os.mkdir(dirname)

counter = 0
cpptraj_script = open('{dirname}/cpptraj.sh'.format(dirname = dirname), 'w')  
for pathway in pathways:
    counter += 1
    cpptraj_in_filename = '{rundir}_{target_bin:05d}_{counter:05d}.cpptraj_in'.format(
                                                                                rundir = rundir,
                                                                                counter = counter,
                                                                                target_bin = target_bin )
    
    cpptraj_in = open('{dirname}/{cpptrajin}'.format(dirname   = dirname,
                                                     cpptrajin = cpptraj_in_filename)
                                                     ,'w')
    for frame in pathway.frame_positions:
        # frame index in cpptraj starts from 1 ( = hdWE segment index + 1) 
        cpptraj_in.write( 'trajin ../{rundir}/{it:05d}.nc {frame} {frame}\n'.format(rundir = rundir, 
                                                                                it     = frame[0],
                                                                                frame  = frame[1] + 1) )
        
    cpptraj_in.write('trajout {rundir}_{target_bin:05d}_{counter:05d}_pathway.nc netcdf'.format(
                                                                                dirname = dirname,
                                                                                rundir = rundir,
                                                                                counter = counter,
                                                                                target_bin = target_bin ))
    cpptraj_in.close() 



    cpptraj_script.write('cpptraj -p $1 -i {cpptraj_in_filename}\n'.format(cpptraj_in_filename = cpptraj_in_filename))
cpptraj_script.close()

# Calculate iteration of divergence matrix
l = len(pathways)
m = numpy.zeros((l,l), dtype = int)
for i in range(l):
    for j in range(l):
        m[i,j] = calcIterationOfDivergence(pathways[i], pathways[j])
numpy.savetxt('{}iteration_of_divergence.dat'.format(dirname), m, fmt = '%05d', 
              header = 'Iteration at which pathway x diverges from pathway y')
for i in range(l):
    m[i,i] = 99999

first_iteration_of_divergence = numpy.min(m)
least_sharing_pathways        = numpy.unravel_index(numpy.argmin(m), m.shape)

for i in range(len(pathways[0].frames)):
    print last_iteration.bins[pathways[0].frames[i][1]].coordinate_ids[0]

sys.stdout.write('\n  Completed.\n'.format(dirname))
sys.stdout.write('   - pathways {p0} and {p1} share least common history, diverging in iteration {it:05d}\n'.format(
                                                            p0 = least_sharing_pathways[0],
                                                            p1 = least_sharing_pathways[1],
                                                            it = first_iteration_of_divergence))
sys.stdout.write('   - shared history information written to {}iteration_of_divergence.dat\n'.format(dirname))
sys.stdout.write('   - cpptraj infiles written to {}*.cpptrajin\n'.format(dirname))
sys.stdout.write('   - cpptraj script written to {}cpptraj.sh\n'.format(dirname))
sys.stdout.write('     (for processing with cpptraj a stripped topology file might be necessary)\n')
sys.stdout.flush()



