import os
import numpy
#~ import multiprocessing as mp
import threading
import sys
from datetime import datetime


class MD_module():
    
    # work_dir                    = ''
    # debug                       = False
    # configuration_file          = ''   path to the amer configuration file that contains the 
    #                                    following entries:
    # amber_topology_path         = ''   path to the amber topology file
    # amber_infile_path           = ''   path to the amber in file
    # amber_coordinate_mask       = ''   for example for RMSD calculation. Example ':1-3@CA'
    # amber_binary                = ''   sander, pmemd, pmemd.cuda
    # parallelization_mode        = ''   serial, parallel, custom_cow
    # number of threads           = ''   number of threads in parallel mode
    
    def __init__(self, work_dir, configuration_file_path, debug):
        """Initializes the MD module and Reads the amber configuration file.
        """
        self.work_dir                 = work_dir
        self.debug                    = debug
        self.configuration_file_path  = configuration_file_path
        
        #read in configuration file
        configuration_file = open(self.configuration_file_path,'r')
        lines              = configuration_file.readlines()
        for line in lines:
            if line.split('=')[0].strip() == 'topology_path':
                self.amber_topology_path   = line.split('=')[1].strip()
            if line.split('=')[0].strip() == 'infile_path':
                self.amber_infile_path   = line.split('=')[1].strip()
            if line.split('=')[0].strip() == 'coordinate_mask':
                self.amber_coordinate_mask   = line.split('=')[1].strip()
            if line.split('=')[0].strip() == 'binary':
                self.amber_binary   = line.split('=')[1].strip()
            if line.split('=')[0].strip() == 'parallelization_mode':
                self.parallelization_mode   = line.split('=')[1].strip()
            if line.split('=')[0].strip() == 'number_of_threads':
                self.number_of_threads   = int(line.split('=')[1].strip())
        configuration_file.close()
    
    def RunMDs(self, iteration):
        """Propagates the trajectories corresponding to an iteration using amber."""
        
        def RunSegmentMD(segment, MD_run_count):
            """Function that runs one single segment MD."""
            command_line = AmberCommandLineString(segment)
            #Command line for debugging
            if self.debug==True:
                    os.system('echo ' + command_line + \
                              ' >> ' + self.work_dir + 'debug/amber_command_lines.log')
            #Log and Run MD
            logfile = open(self.work_dir + 'log/' + iteration.getNameString() + '.MD_log','w')
            logfile.write(MdLogString(segment, status = 0 ))
            sys.stdout.write(writeMdStatus(segment, MD_run_count))
            sys.stdout.flush()
            os.system(command_line)
            logfile.write(MdLogString(segment, status = 1 ))
            logfile.close()
            if self.debug==False:
                RemoveMdOutput(segment)

        def AmberCommandLineString(segment):
            """Returns the command line for an amber run corresponding to the 
            given indices and binary.
            """
            amber_start_coords_path =  self.work_dir + 'run/' + segment.getParentNameString() + '.rst7'
            amber_outfile_path      =  self.work_dir + 'run/' + segment.getNameString()       + '.out'
            amber_trajectory_path   =  self.work_dir + 'run/' + segment.getNameString()       + '.nc'
            amber_end_coords_path   =  self.work_dir + 'run/' + segment.getNameString()       + '.rst7'
        
            amber_command_line      =   self.amber_binary + ' -O' + \
                                        ' -p ' + self.amber_topology_path + \
                                        ' -i ' + self.amber_infile_path + \
                                        ' -c ' + amber_start_coords_path + \
                                        ' -o ' + amber_outfile_path + \
                                        ' -x ' + amber_trajectory_path + \
                                        ' -r ' + amber_end_coords_path
                                        
            return amber_command_line
            
        def RemoveMdOutput(segment):
            """Removes unnecessary MD output files.
            """
            amber_outfile_path      =  self.work_dir + 'run/' + segment.getNameString()       + '.out'
            amber_trajectory_path   =  self.work_dir + 'run/' + segment.getNameString()       + '.nc'
            os.remove(amber_outfile_path)
            os.remove(amber_trajectory_path)
 
            
            

        def MdLogString(segment,status):
            """Returns a string containing system time for MD run logging."""
            if status==0:
                string = 'MD run ' + segment.getNameString() + ' start: ' + str(datetime.now()) + '\n'    
            else:
                string = 'MD run ' + segment.getNameString() + ' end:   ' + str(datetime.now()) + '\n'     
            return string 
            
        def writeMdStatus(segment, MD_run_count):
            """Writes the actual WE run status in a string."""
            number_MD_runs = iteration.getNumberOfSegments()
            string = '\033[1mhdWE Status:\033[0m ' + 'Iteration ' + iteration.getNameString() + \
                     ' Number of bins ' + str(iteration.getNumberOfBins()) + \
                     ' Segment ' + str(MD_run_count).zfill(5) + '/' + str(number_MD_runs).zfill(5) + '\r'
            return string

        
        #Serial Run
        if self.parallelization_mode=='serial':
            MD_run_count = 0
            for bin_loop in iteration:
                for segment_loop in bin_loop:
                    MD_run_count = MD_run_count + 1
                    RunSegmentMD(segment_loop, MD_run_count)
                    
        #Parallel Run   
        if self.parallelization_mode=='parallel':
            MD_run_count = 0
            parallel_jobs = []
            for mBin in iteration:
                for mSegment in mBin:
                    MD_run_count = MD_run_count + 1
                    thread = threading.Thread(target=RunSegmentMD, args=(mSegment, MD_run_count, ))
                    parallel_jobs.append(thread)
                    # TODO 
                    if len(parallel_jobs) >= self.number_of_threads:
                        for job in parallel_jobs:
                            job.start()
                        # Wait until threads are finished
                        for job in parallel_jobs:
                            job.join()
                        # Reset the job list to fill it with next bunch of work
                        parallel_jobs = []
            for job in parallel_jobs:
                job.start()
            # Wait until threads are finished
            for job in parallel_jobs:
                job.join()
            parallel_jobs = []
    
        
    def CalculateCoordinate(self, segment, bins):
        """Calculates the coordinate values of a segment with respect to all
        existing bin reference coordinates and returns them in a numpy array.
        The array entries are in the same order as the bins.
        """

        #Write the cpptraj infile
        segment_name_string=segment.getNameString() 
        
        cpptraj_infile_path=self.work_dir + segment_name_string + '.cpptraj_in'
        cpptraj_infile=open(cpptraj_infile_path,'w')
        cpptraj_infile.write('trajin ' + self.work_dir + 'run/' + segment_name_string + '.rst7' + '\n')     
        cpptraj_output_path = self.work_dir + segment_name_string + '.cpptraj_output'
        for bin_loop in range(0,len(bins)):
            reference_bin_name_string         = self.work_dir + 'run/' + bins[bin_loop].getReferenceNameString() + '.rst7 '
            cpptraj_reference_id_name_string  = '[reference_id_' + str(bin_loop).zfill(5) + ']'

            cpptraj_infile.write('reference ' + reference_bin_name_string + cpptraj_reference_id_name_string + '\n')
            cpptraj_infile.write('rms ' + 'bin_' + str(bin_loop).zfill(5) + ' ' + \
                                 self.amber_coordinate_mask + ' out ' +  cpptraj_output_path + \
                                 ' ref ' + cpptraj_reference_id_name_string + '\n')
        cpptraj_infile.close()


        #Run cpptraj
        cpptraj_execute_string =' -p ' + self.amber_topology_path + \
                                ' -i ' + cpptraj_infile_path + \
                                ' >> ' + self.work_dir + 'debug/cpptraj.log' 
        os.system('cpptraj ' + cpptraj_execute_string )
      
        #Load cpptraj output as numpy array
        try:
            coordinates = numpy.loadtxt(cpptraj_output_path) 
            #Delete the first entry which refers to the frame index
            coordinates = numpy.delete(coordinates, 0)
        except:
            #TODO What should happen then?
            print('amber_module error: cpptraj output ' + cpptraj_output_path + ' can not be found or loaded.')
        
        #Check if cpptraj output is correct
        if not (len(coordinates)==len(bins)):
            #TODO What should happen then?
            print('amber_module error: cpptraj output ' + cpptraj_output_path + ' does not have the correct number of entries.' )

        #Remove temporary files
        if (self.debug==False):
            os.remove(cpptraj_infile_path)
            os.remove(cpptraj_output_path)
            os.remove(self.work_dir + 'debug/cpptraj.log')
        
        return coordinates
        
    def ana_calculatePMF_getCoordinate(self, segment, cpptraj_lines):
        """
        Calculates the value of a coordinate corresponding to a segment and defined
        in cpptraj_line via cpptraj.
        """
        
        #Write the cpptraj infile
        segment_name_string=segment.getNameString() 
        cpptraj_infile_path=self.work_dir + segment_name_string + '.ana_calculatePMF_cpptraj_in'
        cpptraj_output_path = self.work_dir + segment_name_string + '.ana_calculatePMF_cpptraj_output'
        cpptraj_infile=open(cpptraj_infile_path,'w')
        cpptraj_infile.write('trajin ' + self.work_dir + 'run/' + segment_name_string + '.rst7' + '\n')
        cpptraj_infile.writelines(cpptraj_lines + ' out ' + cpptraj_output_path )
        cpptraj_infile.close()
        
        #Execute cpptraj
        cpptraj_execute_string =' -p ' + self.amber_topology_path + \
                                ' -i ' + cpptraj_infile_path + \
                                ' >> ' + self.work_dir + 'debug/ana_calculatePMF_cpptraj.log' 
        os.system('cpptraj ' + cpptraj_execute_string )        
        
        #Load cpptraj output as numpy array
        try:
            coordinates = numpy.loadtxt(cpptraj_output_path) 
        except:
            #TODO What should happen then?
            print('amber_module error: cpptraj output ' + cpptraj_output_path + ' can not be found or loaded.')
            
        #Remove temporary files
        if (self.debug==False):
            os.remove(cpptraj_infile_path)
            os.remove(cpptraj_output_path)
            os.remove(self.work_dir + 'debug/ana_calculatePMF_cpptraj.log' )
            
        coordinate_value = coordinates[1]
        
        return coordinate_value
