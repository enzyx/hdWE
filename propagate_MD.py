import os

def propagate_AMBER(directory,IT,BIN,SEG,dir_infile,dir_topology,remove_nonobligatory_md_output):
    """Propagates the trajectory corresponding to given iteration, bin and segment
    with the pmemd module of AMBER."""

       
    string_trajectory_index =   'IT' + str(IT).zfill(4) + \
                                '_BIN' + str(BIN).zfill(4) + \
                                '_SEG' + str(SEG).zfill(4)
                                
    dir_start_coords        =   directory + '/' + string_trajectory_index+'.rst7'
    dir_end_coords          =   directory + '/' + string_trajectory_index+'_end.rst7'
    dir_trajectory_coords   =   directory + '/' + string_trajectory_index+'.nc'
    dir_outfile             =   directory + '/' + string_trajectory_index+'.out'
    
    exec_amber_string       =   'pmemd -O' + \
                                ' -p ' + dir_topology + \
                                ' -i ' + dir_infile + \
                                ' -c ' + dir_start_coords + \
                                ' -o ' + dir_outfile + \
                                ' -x ' + dir_trajectory_coords + \
                                ' -r ' + dir_end_coords
    
    os.system('echo ' + exec_amber_string + ' >> ' + directory + '/MD_AMBER_command_line.log')
    os.system(exec_amber_string)
    
    if remove_nonobligatory_md_output==True:
        os.remove(dir_trajectory_coords)
        os.remove(dir_outfile)    
    
                            
    
    
    
