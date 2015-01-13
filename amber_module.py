import os
import numpy

class MD_module():
    
    work_dir                    = ''
    debug                       = False
    configuration_file          = ''  # path to the amer configuration file that contains the 
                                      # following entries:
    amber_topology_path         = ''  # path to the amber topology file
    amber_infile_path           = ''  # path to the amber in file
    amber_coordinate_mask       = ''  # for example for RMSD calculation. Example ':1-3@CA'
    amber_binary                = ''  # sander, pmemd, pmemd.cuda
    parallelization_mode        = ''  # serial_cpu, parallel_cpu, cuda, custom_cow, add a gusto
    
    def __init__(self, configuration_file_path, work_dir, debug):
        """Initializes the MD module and Reads the amber configuration file.
        """
        self.work_dir                 = work_dir
        self.debug                    = debug
        self.configuration_file_path  = configuration_file_path
        
        #read in configuration file
        #TODO
        configuration_file         =open(self.configuration_file_path,'r')
        self.amber_topology_path   = configuration_file.readline().strip()
        self.amber_infile_path     = configuration_file.readline().strip()
        self.amber_coordinate_mask = configuration_file.readline().strip()
        self.amber_binary          = configuration_file.readline().strip()
        self.parallelization_mode  = configuration_file.readline().strip()
        configuration_file.close()
    
    def RunMDs(self, iteration):
        """Propagates the trajectories corresponding to an iteration using amber.
        """

        def AmberCommandLineString(segment):
            """Returns the command line for an amber run corresponding to the 
            given indices.
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
            
        if self.parallelization_mode=='serial_cpu':
            for bin_loop in iteration.bins:
                for segment_loop in bin_loop.segments:
                    if self.debug==True:
                         command_line = AmberCommandLineString(segment_loop)
                         os.system('echo ' + command_line + \
                                   ' >> ' + self.work_dir + 'debug/MD_AMBER_command_lines.log')
                    os.system(command_line)
    

        
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
            print('amber_module error: cpptraj output ' + cpptraj_output_path + ' can not be found or loaded.')
        
        #Check if cpptraj output is correct
        if not (len(coordinates)==len(bins)):
            print('amber_module error: cpptraj output ' + cpptraj_output_path + ' does not have the correct number of entries.' )

        #Remove temporary files
        if (self.debug==False):
            os.remove(cpptraj_infile_path)
            os.remove(cpptraj_output_path)
            os.remove(self.work_dir + 'debug/cpptraj.log')
        
        return coordinates    
        
