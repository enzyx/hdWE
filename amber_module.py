from __future__ import print_function
import os
import ConfigParser
import numpy
import threading
import sys
from datetime import datetime
from thread_container import ThreadContainer
import pickle

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
       
        # check topology and infile:
        if not os.path.isfile(self.amber_topology_path):
            raise Exception("No topology found at given path.")
        if not os.path.isfile(self.amber_infile_path):
            raise Exception("No infile found at given path.")        
    
    def setIteration(self, iteration):
        self.iteration = iteration
    
    def RunSegmentMD(self, segment, MD_run_count, MD_skip_count):
        """Function that runs one single segment MD."""
        command_line = self.AmberCommandLineString(segment)
        #Command line for debugging
        if self.debug==True:
                os.system('echo ' + command_line + \
                          ' >> ' + self.workdir + 'debug/amber_command_lines.log')
        #Log and Run MD
        logfile = open(self.workdir + 'log/' + self.iteration.getNameString() + '.MD_log','a')
        logfile.write(self.MdLogString(segment, status = 0 ))
        sys.stdout.write(self.writeMdStatus(segment, MD_run_count, MD_skip_count))
        sys.stdout.flush()
        os.system(command_line)
        logfile.write(self.MdLogString(segment, status = 1 ))
        logfile.close()
        if not self.debug:
            self.RemoveMdOutput(segment)

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
    
    def AmberCommandLineString(self, segment):
        """Returns the command line for an amber run corresponding to the 
        given indices and binary.
        """
        amber_start_coords_path = self.workdir + 'run/' + segment.getParentNameString() + '.rst7'
        amber_outfile_path      = self.workdir + 'run/' + segment.getNameString()       + '.out'
        amber_trajectory_path   = self.workdir + 'run/' + segment.getNameString()       + '.nc'
        amber_end_coords_path   = self.workdir + 'run/' + segment.getNameString()       + '.rst7'
    
        amber_command_line      = self.amber_binary + ' -O' + \
                                  ' -p ' + self.amber_topology_path + \
                                  ' -i ' + self.amber_infile_path + \
                                  ' -c ' + amber_start_coords_path + \
                                  ' -o ' + amber_outfile_path + \
                                  ' -x ' + amber_trajectory_path + \
                                  ' -r ' + amber_end_coords_path
                                    
        return amber_command_line
    
    def SkipCommandLineString(self, segment):
        """Returns the command line for linking segment restart files of skipped bins to next iteration.
        """
        amber_start_coords_path = self.workdir + 'run/' + segment.getParentNameString() + '.rst7'
        amber_end_coords_path   = self.workdir + 'run/' + segment.getNameString()       + '.rst7'
        
        skip_command_line       = 'ln -s' + \
                                  '  ' + amber_start_coords_path + \
                                  '  ' + amber_end_coords_path
        return skip_command_line
    
        
    def RemoveMdOutput(self, segment):
        """Removes unnecessary MD output files.
        """
        amber_outfile_path      =  self.workdir + 'run/' + segment.getNameString()       + '.out'
        amber_trajectory_path   =  self.workdir + 'run/' + segment.getNameString()       + '.nc'
        os.remove(amber_outfile_path)
        os.remove(amber_trajectory_path)
    
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
            for bin_loop in iteration:
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
            for bin_loop in iteration:
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
            iteration_dump_filename = self.workdir + 'run/' + 'iteration.dump.tmp'
            iteration_dump_file = open(iteration_dump_filename, 'w')
            pickle.dump(iteration, iteration_dump_file, protocol=2)
            iteration_dump_file.close()
            # Call myself with MPI...
            if self.debug:
                print("Switching to MPI for MD now...")
            os.system(self.mpirun + ' ' + 
                      sys.executable + ' ' +
                      os.path.realpath(__file__) + ' ' + 
                      self.configfile + ' ' + 
                      str(self.debug) + ' ' + 
                      iteration_dump_filename)

        
    def calcRmsdToBins(self, segment, bins):
        """
        Calculates the coordinate values of a segment with respect to all
        existing bin reference coordinates and returns them in a numpy array.
        The array entries are in the same order as the bins.
        """

        #Write the cpptraj infile
        segment_name_string = segment.getNameString() 
        
        cpptraj_infile_path = "{workdir}/{segment}.cpptraj_in".format(workdir=self.workdir, 
                                                                      segment=segment_name_string)
        cpptraj_infile      = open(cpptraj_infile_path, 'w')
        cpptraj_infile.write('trajin {workdir}/run/{segment}.rst7\n'.format(workdir=self.workdir, 
                                                                            segment=segment_name_string))  
        cpptraj_output_path = "{workdir}/{segment}.cpptraj_output".format(workdir=self.workdir, 
                                                                          segment=segment_name_string)   
        
        for bin_loop in range(0,len(bins)):
            reference_bin_name_string         = "{workdir}/run/{segment}.rst7".format(workdir=self.workdir, 
                                                                                      segment=bins[bin_loop].getReferenceNameString())
            cpptraj_reference_id_name_string  = '[reference_id_{0:05d}]'.format(bin_loop)

            cpptraj_infile.write('reference {0} {1}\n'.format(reference_bin_name_string, cpptraj_reference_id_name_string))
            cpptraj_infile.write('rms in_{0:05d} {1} out {2} ref {3}\n'.format(bin_loop,
                                 self.amber_coordinate_mask, cpptraj_output_path,
                                 cpptraj_reference_id_name_string))
        cpptraj_infile.close()

        #Run cpptraj
        cpptraj_execute_string = ' -p {top} -i {inpath} >> {workdir}/log/cpptraj.log'.format(
                                                            top=self.amber_topology_path, 
                                                            inpath=cpptraj_infile_path,
                                                            workdir=self.workdir)
        os.system('cpptraj {execute}'.format(execute=cpptraj_execute_string))
      
        #Load cpptraj output as numpy array
        try:
            coordinates = numpy.loadtxt(cpptraj_output_path) 
            #Delete the first entry which refers to the frame index
            coordinates = numpy.delete(coordinates, 0)
        except:
            #TODO What should happen then?
            print('amber_module error: cpptraj output {0} can not '\
                  'be found or loaded.'.format(cpptraj_output_path))
        
        #Check if cpptraj output is correct
        if not (len(coordinates)==len(bins)):
            #TODO What should happen then?
            print("amber_module error: cpptraj output {0} does not " \
                  "have the correct number of entries.".format(cpptraj_output_path))

        #Remove temporary files
        if not self.debug:
            os.remove(cpptraj_infile_path)
            os.remove(cpptraj_output_path)
        
        return coordinates
    
    def ana_calcCoordinateOfSegment(self, segment, cpptraj_lines):
        """
        Calculates the value of a coordinate corresponding to a segment and defined
        in cpptraj_line via cpptraj.
        """
        
        #Write the cpptraj infile
        segment_name_string=segment.getNameString() 
        cpptraj_infile_path=self.workdir + segment_name_string + '.ana_calculatePMF_cpptraj_in'
        cpptraj_output_path = self.workdir + segment_name_string + '.ana_calculatePMF_cpptraj_output'
        cpptraj_infile=open(cpptraj_infile_path,'w')
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

# Call only in MPI mode
if __name__ == "__main__":
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # Read in the call arguments
    configfile               = sys.argv[1]
    debug                    = bool(str(sys.argv[2]) == "True")
    iteration_dump_file_path = sys.argv[3]
    # Initialize MD module
    md_module = MD_module(configfile = configfile, debug = debug)
    iteration_dump_file = open(iteration_dump_file_path, "r")
    md_module.setIteration(pickle.load(iteration_dump_file))
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
            os.remove(iteration_dump_file_path)
