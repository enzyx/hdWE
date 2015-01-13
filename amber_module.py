import os
import numpy

class MD_module():
    
    work_dir                    = ''
    configuration_file          = ''  # path to the amer configuration file that contains the 
                                      # following entries:
    amber_topology_path         = ''  # path to the amber topology file
    amber_infile_path           = ''  # path to the amber in file
    amber_coordinate_mask       = ''  # for example for RMSD calculation. Example ':1-3@CA'
    amber_binary                = ''  # sander, pmemd, pmemd.cuda
    parallelization_mode        = ''  # serial_cpu, parallel_cpu, cuda, custom_cow, add a gusto
    
    def __init__(self, configuration_file_path, work_dir):
        """Initializes the MD module and Reads the amber configuration file.
        """
        self.work_dir                 = work_dir
        self.configuration_file_path   = configuration_file_path
        
        #read in configuration file
        #TODO
    
    def RunMDs(self, work_dir,MD_mode,debug,trajectory_index,trajectory_parent_index):
        """Propagates the trajectory corresponding to given index with AMBER."""

        sub_dir = 'propagation/'
        
        if MD_mode=='pmemd':
            binary='pmemd'
        elif MD_mode=='sander':
            binary='sander'
        elif MD_mode=='pmemd.cuda':
            binary='pmemd.cuda'
        else:
            print('Error: MD mode {mode:s} not known.'.format(mode=MD_mode))

           
        string_trajectory_index =   str(IT).zfill(5) + \
                                    str(BIN).zfill(5) + \
                                    str(SEG).zfill(5)
                                    
        dir_topology            =   work_dir + 'MD_infiles/system_amber.topology'
        dir_infile              =   work_dir + 'MD_infiles/system_amber.mdin'
        dir_start_coords        =   work_dir + sub_dir + string_trajectory_parent_index+'.rst7'
        dir_end_coords          =   work_dir + sub_dir + string_trajectory_index+'.rst7'
        dir_trajectory_coords   =   work_dir + sub_dir + string_trajectory_index+'.nc'
        dir_outfile             =   work_dir + sub_dir + string_trajectory_index+'.out'
        
        exec_amber_string       =   binary + ' -O' + \
                                    ' -p ' + dir_topology + \
                                    ' -i ' + dir_infile + \
                                    ' -c ' + dir_start_coords + \
                                    ' -o ' + dir_outfile + \
                                    ' -x ' + dir_trajectory_coords + \
                                    ' -r ' + dir_end_coords
        
        os.system('echo ' + exec_amber_string + ' >> ' + work_dir + 'debug/MD_AMBER_command_lines.log')
        os.system(exec_amber_string)
    
        if debug==False:
            os.remove(dir_trajectory_coords)
            os.remove(dir_outfile)    
        
    def CalculateCoordinate(self, segment, bins):
        """Calculates the coordinate values of a segment with respect to all
        existing bin reference coordinates and returns them in an array.
        The array entries are in the same order as the bins.
        """
        segment_name_string=segment.getNameString() 

        #Write the cpptraj infile
        cpptraj_infile=open(self.work_dir + 'cpptraj_' + segment_name_string + '.infile','w')
        cpptraj_infile.write('trajin ' + segment_name_string + '.rst7')       
        for bin_loop in range(0,len(bins)):

            cpptraj_infile.write('reference ' + self.work_dir + 'run/' + bins[bin_loop].getNameString() + 
                                '.rst7 reference_id_' + str(bin_loop).zfill(5) )
            cpptraj_infile.write('rms  + reference_id_' + bins[bin_loop].getNameString() + ' ' + 
                                self.coordinate_mask + ' out ' +  self._work_dir + 
                                'cpptraj_'+segment_name_string + '.output' )
        
        cpptraj_infile.close()
        #Run cpptraj
        cppraj_execute_string = ' -p ' + self.amber_topology_path + \
                                ' -i ' + self.amber_topology_path + \
                                ' -log ' + self.work_dir + '/debug/cpptraj_' + segment_name_string + '.log'
        os.system('cpptraj ' + cpptraj_execute_string)
      
        #Load cpptraj output as numpy array
        try:
            coordinates = numpy.loadtxt(self.work_dir + 'cpptraj_' + segment_name_string + '.output') 
        except:
            print('cpptraj output can not be found or loaded.')
            #TODO what should happen then?
        
        #Check if cpptraj output is correct
        if not (len(coordinates)==len(bins)):
            #TODO
        
        return coordinates    
        
