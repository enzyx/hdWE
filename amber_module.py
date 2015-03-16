from __future__ import print_function
import os, shutil
import ConfigParser
import numpy
import threading
import sys
from datetime import datetime
from thread_container import ThreadContainer
import pickle
import tempfile
import uuid

class MD_module():
    
    # workdir                    = ''
    # debug                       = False
    # configuration_file          = ''   path to the amer configuration file that contains the 
    #                                    following entries:
    # amber_topology_path         = ''   path to the amber topology file
    # amber_infile_path           = ''   path to the amber in file
    # amber_coordinate_mask       = ''   for example for RMSD calculation. Example ':1-3@CA'
    # amber_binary                = ''   sander, pmemd, pmemd.cuda
    # parallelization_mode        = ''   serial, parallel, custom_cow
    # number of threads           = ''   number of threads in parallel mode
    
    def __init__(self, configfile, debug):
        """Initializes the MD module and Reads the amber configuration file.
        """
        #read in configuration file
        self.configfile = configfile
        config = ConfigParser.ConfigParser()
        config.read(configfile)
        self.workdir               = config.get('hdWE','workdir')
        self.amber_topology_path   = self.workdir + config.get('amber','topology-path')
        self.amber_infile_path     = self.workdir + config.get('amber','infile-path')
        self.amber_coordinate_mask = config.get('amber','coordinate-mask')
        self.amber_binary          = config.get('amber','binary')
        self.parallelization_mode  = config.get('amber','parallelization-mode')
        # Only search mpirun binary in config file if MPI switched on
        if self.parallelization_mode == 'mpi':
            self.mpirun                = config.get('amber', 'mpirun')
        self.number_of_threads     = int(config.get('amber','number-of-threads'))
        self.debug                 = debug
        # local link to iteration
        self.iteration             = None
        self.ITERATION_DUMP_FNAME  = self.workdir + '/run/iteration.dump'
        self.RMSD_MATRIX_DUMP_FNAME= self.workdir + '/run/rmsd_matrix.dump'
                       
        # check topology and infile:
        if not os.path.isfile(self.amber_topology_path):
            raise Exception("No topology found at given path.")
        if not os.path.isfile(self.amber_infile_path):
            raise Exception("No infile found at given path.")        
    
    def setIteration(self, iteration):
        self.iteration = iteration
    
    def RunSegmentMD(self, segment, MD_run_count, MD_skip_count):
        """
        Function that runs one single segment MD.
        """
        amber_start_coords_path = "{0}/run/{1}.rst7".format(self.workdir, segment.getParentNameString()) 
        amber_end_coords_path   = "{0}/run/{1}.rst7".format(self.workdir, segment.getNameString()) 
        
        if self.debug:
            amber_info_path         = "{0}/run/{1}.inf".format(self.workdir, segment.getNameString())
            amber_outfile_path      = "{0}/run/{1}.out".format(self.workdir, segment.getNameString())
            amber_trajectory_path   = "{0}/run/{1}.nc".format(self.workdir, segment.getNameString())
        else:
            amber_info_path         = '/dev/null'
            amber_trajectory_path   = '/dev/null'
            #amber_outfile_path      = '/dev/null'
            #TODO: Strange bug on REX allows only info and traj to be /dev/null
            #amber_outfile    = tempfile.NamedTemporaryFile(suffix=segment.getNameString(), prefix=".out", delete=True)
            #amber_outfile_path    = amber_outfile.name
            amber_outfile_path     = '/tmp/{0}_{1}.out'.format(segment.getNameString(), uuid.uuid1()) 
        
        command_line = self.amber_binary + ' -O' + \
                                  ' -p '   + self.amber_topology_path + \
                                  ' -i '   + self.amber_infile_path + \
                                  ' -c '   + amber_start_coords_path + \
                                  ' -o '   + amber_outfile_path + \
                                  ' -x '   + amber_trajectory_path + \
                                  ' -inf ' + amber_info_path + \
                                  ' -r '   + amber_end_coords_path
                                    
        

        #Command line for debugging
        if self.debug:
            os.system('echo {0} >> {1}/debug/amber_command_lines.log'.format(command_line, self.workdir))
                      
        #Log and Run MD
        if self.debug:
            logfile = open("{0}/log/{1}.MD_log".format(self.workdir, self.iteration.getNameString()),'a')
            logfile.write(self.MdLogString(segment, status = 0))
        sys.stdout.write(self.writeMdStatus(segment, MD_run_count, MD_skip_count))
        sys.stdout.flush()
        # Run MD
        os.system(command_line)
        if self.debug:
            logfile.write(self.MdLogString(segment, status = 1 ))
            logfile.close()
        
        # Remove out file
        if not self.debug:
            try: os.remove(amber_outfile_path)
            except OSError: pass
        
    def SkipSegmentMD(self, segment, MD_run_count, MD_skip_count):
        """Function that runs one single segment MD."""
        command_line = self.SkipCommandLineString(segment)
        #Command line for debugging
        if self.debug:
                os.system('echo ' + command_line + \
                          ' >> ' + self.workdir + 'debug/amber_command_lines.log')
        #Log and Run MD
        logfile = open(self.workdir + 'log/' + self.iteration.getNameString() + '.MD_log','a')
        logfile.write(self.MdLogString(segment, status = 2 ))
        sys.stdout.write(self.writeMdStatus(segment, MD_run_count, MD_skip_count))
        sys.stdout.flush()
        os.system(command_line)
        logfile.close()
    
    def SkipCommandLineString(self, segment):
        """
        Returns the command line for linking segment restart files of skipped bins to next iteration.
        """
        amber_start_coords_path = self.workdir + 'run/' + segment.getParentNameString() + '.rst7'
        amber_end_coords_path   = self.workdir + 'run/' + segment.getNameString()       + '.rst7'
        
        skip_command_line       = 'ln {0} {1}'.format(
                                                  amber_start_coords_path,
                                                  amber_end_coords_path)
        return skip_command_line
        
    def MdLogString(self, segment,status):
        """Returns a string containing system time for MD run logging."""
        if status==0:
            string = 'MD run ' + segment.getNameString() + ' start: ' + str(datetime.now()) + '\n'    
        elif status==1:
            string = 'MD run ' + segment.getNameString() + ' end:   ' + str(datetime.now()) + '\n'     
        else:
            string = 'MD run ' + segment.getNameString() + ' skipped: ' + str(datetime.now()) + '\n'    
        return string 
        
    def writeMdStatus(self, segment, MD_run_count, MD_skip_count):
        """Writes the actual WE run status in a string."""
        number_MD_runs = self.iteration.getNumberOfSegments()
        string = '\033[1mhdWE Status:\033[0m ' + 'Iteration ' + self.iteration.getNameString() + \
                 ' Number of bins ' + str(self.iteration.getNumberOfBins()) + \
                 ' Segment ' + str(MD_run_count).zfill(5) + '/' + str(number_MD_runs).zfill(5) + \
                 ' Skipped segments: ' + str(MD_skip_count).zfill(6) + '\r'
        return string
    
    def RunMDs(self, iteration):
        """Propagates the trajectories corresponding to an iteration using amber."""
        self.iteration = iteration
        #Serial Run
        if self.parallelization_mode=='serial':
            MD_run_count  = 0
            MD_skip_count = 0
            for bin_loop in self.iteration:
                if bin_loop.isConverged() == False:
                    for segment_loop in bin_loop:
                        MD_run_count  += 1
                        self.RunSegmentMD(segment_loop, MD_run_count, MD_skip_count)
                else:
                    for segment_loop in bin_loop:
                        MD_skip_count += 1  
                        MD_run_count  += 1
                        self.SkipSegmentMD(segment_loop, MD_run_count, MD_skip_count)                    
                    
        #Thread Parallel Run   
        if self.parallelization_mode=='thread':
            MD_run_count  = 0
            MD_skip_count = 0
            thread_container = ThreadContainer()
            for bin_loop in self.iteration:
                if bin_loop.isConverged() == False:
                    for segment_loop in bin_loop:
                        MD_run_count  += 1
                        thread_container.appendJob(threading.Thread(target=self.RunSegmentMD, 
                                                                    args=(segment_loop, MD_run_count, MD_skip_count, )))
                        if thread_container.getNumberOfJobs() >= self.number_of_threads:
                            thread_container.runJobs()
                else:
                    for segment_loop in bin_loop:
                        MD_skip_count += 1
                        MD_run_count  += 1
                        self.SkipSegmentMD(segment_loop, MD_run_count, MD_skip_count)     
                        
            # Finish jobs in queue
            thread_container.runJobs()
        
        #MPI parallel run
        if self.parallelization_mode == 'mpi':
            # dump the iteration object to a file
            self.dumpIterationToFile()
            # Call myself with MPI...
            if self.debug:
                print("Switching to MPI for MD now...")
            os.system("{0} {1} {2} {3} {4} {5} ".format(self.mpirun, 
                                  sys.executable,
                                  os.path.realpath(__file__),
                                  "MPIMD",
                                  self.configfile,
                                  str(self.debug)))

        
    def calcRmsdToBins(self, segment, bins):
        """
        Calculates the coordinate values of a segment with respect to all
        existing bin reference coordinates and returns them in a numpy array.
        The array entries are in the same order as the bins.
        """

        #Write the cpptraj infile
        segment_name_string = segment.getNameString() 
        UUID = uuid.uuid1()
        
        if self.debug:
            cpptraj_infile_path = "{workdir}/{segment}.cpptraj_in".format(workdir=self.workdir, 
                                                                      segment=segment_name_string)
        else:
            cpptraj_infile_path = "/tmp/{0}_{1}.cpptraj_in".format(segment_name_string, UUID)
        
        cpptraj_infile      = open(cpptraj_infile_path, 'w')
        cpptraj_infile.write('trajin {workdir}/run/{segment}.rst7\n'.format(workdir=self.workdir, 
                                                                            segment=segment_name_string))
        if self.debug:  
            cpptraj_outfile_path = "{workdir}/{segment}.cpptraj_out".format(workdir=self.workdir, 
                                                                      segment=segment_name_string)
        else:
            cpptraj_outfile_path = "/tmp/{0}_{1}.cpptraj_out".format(segment_name_string, UUID)
        
        for bin_loop in range(0,len(bins)):
            reference_bin_name_string         = "{workdir}/run/{segment}.rst7".format(workdir=self.workdir, 
                                                                                      segment=bins[bin_loop].getReferenceNameString())
            cpptraj_reference_id_name_string  = '[reference_id_{0:05d}]'.format(bin_loop)

            cpptraj_infile.write('reference {0} {1}\n'.format(reference_bin_name_string, cpptraj_reference_id_name_string))
            cpptraj_infile.write('rms in_{0:05d} {1} out {2} ref {3}\n'.format(bin_loop,
                                 self.amber_coordinate_mask, cpptraj_outfile_path,
                                 cpptraj_reference_id_name_string))
        # Write changes to file
        cpptraj_infile.close()

        #Run cpptraj
        if self.debug:
            cpptraj_execute_string = ' -p {top} -i {inpath} > {workdir}/log/cpptraj.log'.format(
                                                            top=self.amber_topology_path, 
                                                            inpath=cpptraj_infile_path,
                                                            workdir=self.workdir)
        else:
            cpptraj_execute_string = ' -p {top} -i {inpath} > /dev/null'.format(
                                                            top=self.amber_topology_path, 
                                                            inpath=cpptraj_infile_path)
        os.system('cpptraj {execute}'.format(execute=cpptraj_execute_string))
                
        #Load cpptraj output as numpy array
        try:
            coordinates = numpy.loadtxt(cpptraj_outfile_path) 
            #Delete the first entry which refers to the frame index
            coordinates = numpy.delete(coordinates, 0)
        except:
            #TODO What should happen then?
            print('amber_module error: cpptraj output {0} can not '\
                  'be found or loaded.'.format(cpptraj_outfile_path))
            
        #Check if cpptraj output is correct
        if not (len(coordinates)==len(bins)):
            #TODO What should happen then?
            print("amber_module error: cpptraj output {0} does not " \
                  "have the correct number of entries.".format(cpptraj_outfile_path))

        #Remove temporary files
        if not self.debug: 
            try: 
                os.remove(cpptraj_outfile_path)
                os.remove(cpptraj_infile_path)
            except OSError: pass
        
        return coordinates
    
    def loadIterationFromDumpFile(self):
        """
        Load the local copy of iteration from the 
        default dump file name 
        """
        iteration_dump_file = open(self.ITERATION_DUMP_FNAME, "r")
        self.setIteration(pickle.load(iteration_dump_file))
        iteration_dump_file.close()
    
    def dumpIterationToFile(self):
        """
        Dump the local iteration object to file
        """
        iteration_dump_file = open(self.ITERATION_DUMP_FNAME, 'w')
        pickle.dump(self.iteration, iteration_dump_file, protocol=2)
        iteration_dump_file.close()
    
    def removeIterationDumpFile(self):
        """
        Delete the default iteration dump file
        """
        os.remove(self.ITERATION_DUMP_FNAME)
    
    def dumpRmsdMatrixToFile(self, rmsd_matrix):
        """
        Store the rmsd matrix to dump file
        """        
        rmsd_matrix_dump_file = open(self.RMSD_MATRIX_DUMP_FNAME , 'w')
        rmsd_matrix.dump(rmsd_matrix_dump_file)
        rmsd_matrix_dump_file.close()

    def loadRmsdMatrixFromDumpFile(self):
        """
        Load the rmsd matrix from dump file
        """
        rmsd_matrix_dump_file = open(self.RMSD_MATRIX_DUMP_FNAME , 'r')
        rmsd_matrix = pickle.load(rmsd_matrix_dump_file)
        rmsd_matrix_dump_file.close()
        return rmsd_matrix

    def removeRmsdMatrixDumpFile(self):
        """ Remove rmsd matrix dump file """
        os.remove(self.RMSD_MATRIX_DUMP_FNAME)
    
    def calcRmsdSegmentsToBinsMatrix(self, iteration):
        """
        Returns a NxM matrix with RMSD entries segments x bins
        """
        self.setIteration(iteration)
        if self.parallelization_mode in ["serial", "thread"]:
            rmsd_matrix = numpy.zeros((iteration.getNumberOfSegments(), iteration.getNumberOfBins()))
            # calculate entries
            i = 0
            for bin_loop in self.iteration:
                for segment in bin_loop:
                    coordinates = self.calcRmsdToBins(segment, self.iteration.bins)
                    # fill matrix
                    for j in range(len(coordinates)):
                        rmsd_matrix[i][j] = coordinates[j]
                    i += 1
            return rmsd_matrix
        
        if self.parallelization_mode == "mpi":            
            # dump the iteration object to a file
            self.dumpIterationToFile()
            # Call myself with MPI...
            if self.debug:
                print("Switching to MPI for RMSD calculation.")
            os.system("{0} {1} {2} {3} {4} {5} ".format(self.mpirun, 
                                  sys.executable,
                                  os.path.realpath(__file__),
                                  "MPIRMSD",
                                  self.configfile, 
                                  str(self.debug)))
            rmsd_matrix = self.loadRmsdMatrixFromDumpFile()
            if not self.debug:
                self.removeRmsdMatrixDumpFile()
            return rmsd_matrix
    
    def ana_calcCoordinateOfSegment(self, segment, cpptraj_lines):
        """
        Calculates the value of a coordinate corresponding to a segment and defined
        in cpptraj_line via cpptraj.
        """
        
        #Write the cpptraj infile
        segment_name_string = segment.getNameString() 
        cpptraj_infile_path = self.workdir + segment_name_string + '.ana_calculatePMF_cpptraj_in'
        cpptraj_output_path = self.workdir + segment_name_string + '.ana_calculatePMF_cpptraj_output'
        cpptraj_infile      = open(cpptraj_infile_path,'w')
        cpptraj_infile.write('trajin ' + self.workdir + 'run/' + segment_name_string + '.rst7' + '\n')
        cpptraj_infile.writelines(cpptraj_lines + ' out ' + cpptraj_output_path )
        cpptraj_infile.close()
        
        #Execute cpptraj
        cpptraj_execute_string =' -p ' + self.amber_topology_path + \
                                ' -i ' + cpptraj_infile_path
        cpptraj_execute_string = cpptraj_execute_string + ' >> ' + self.workdir + 'log/ana_calculatePMF_cpptraj.log' 
        os.system('cpptraj ' + cpptraj_execute_string )        
        
        #Load cpptraj output as numpy array
        try:
            coordinates = numpy.loadtxt(cpptraj_output_path) 
        except:
            #TODO What should happen then?
            print('amber_module error: cpptraj output ' + cpptraj_output_path + ' can not be found or loaded.')
            
        #Remove temporary files
        if not self.debug:
            os.remove(cpptraj_infile_path)
            os.remove(cpptraj_output_path)
            
        coordinate_value = coordinates[1]
        
        return coordinate_value


def doMPIMD(configfile, debug):
    """
    This function can run on multiple MPI processes in parallel
    to propagate all segments of iteration via MD 
    """    
    # Read in the call arguments
    configfile               = configfile
    debug                    = bool(debug == "True")
    # Initialize MD module
    md_module = MD_module(configfile = configfile, debug = debug)
    md_module.loadIterationFromDumpFile()
    if debug:
        if rank == 0:
            sys.stderr.write("Number of MPI processes: {0}\n".format(size))
            sys.stderr.flush()
        comm.barrier()
    # Every node works only on segments modulo their rank 
    workcount = 0
    md_skip_count = 0
    for loop_bin in md_module.iteration:
        for loop_segment in loop_bin:
            if not loop_bin.isConverged():
                if workcount % size == rank:
                    # Run MD on this node
                    md_module.RunSegmentMD(loop_segment, workcount, md_skip_count)
            else:
                md_skip_count += 1
                if workcount % size == rank:
                    # Run MD skip
                    md_module.SkipSegmentMD(loop_segment, workcount, md_skip_count)
            workcount += 1
    if rank == 0:
        sys.stderr.write("\n")
        if debug:
            sys.stderr.write("Finishing MPI\n")
        sys.stderr.flush()
        # Remove the iteration dump file
        if not debug:
            md_module.removeIterationDumpFile()

def doMPICalcRmsdMatrix(configfile, debug):
    # Read in the call arguments
    configfile               = configfile
    debug                    = bool(debug == "True")
    # Initialize MD module
    md_module = MD_module(configfile = configfile, debug = debug)
    md_module.loadIterationFromDumpFile()
    iteration = md_module.iteration
    
    if debug and rank == 0:
        sys.stderr.write("Number of MPI processes: {0}\n".format(size))
        sys.stderr.flush()

    # Setup the rmsd matrix
    rmsd_matrix = numpy.zeros((iteration.getNumberOfSegments(), iteration.getNumberOfBins()))
    # calculate entries
    i = 0
    for bin_loop in iteration:
        for segment in bin_loop:
            if i % size == rank:
                coordinates = md_module.calcRmsdToBins(segment, iteration.bins)
                # fill matrix
                for j in range(len(coordinates)):
                    rmsd_matrix[i][j] = coordinates[j]
            i += 1
    # Collect the matrix on the root process (rank == 0)
    rmsd_matrix = comm.reduce(rmsd_matrix, op = MPI.SUM, root=0)
    # Dump matrix to file
    if rank == 0:
        md_module.dumpRmsdMatrixToFile(rmsd_matrix)
        if debug:
            print("RMSD Matrix: ", rmsd_matrix)
            print("Finished MPI RMSD Matrix calculation")

    

# Call only in MPI mode
if __name__ == "__main__":
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    if sys.argv[1] == "MPIMD":
        doMPIMD(configfile=sys.argv[2], 
                debug=sys.argv[3])
    if sys.argv[1] == "MPIRMSD":
        doMPICalcRmsdMatrix(configfile=sys.argv[2],
                            debug=sys.argv[3])
        