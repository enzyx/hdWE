import os

def run_MD(work_dir,MD_mode,debug,trajectory_index,trajectory_parent_index):
    """Propagates the trajectory corresponding to given index with the pmemd module of AMBER."""

    sub_dir = 'propagation/'
    
    if MD_mode=='pmemd':
        binary='pmemd'
    elif MD_mode=='sander':
        binary='sander'
    elif MD_mode=='pmemd.cuda':
        binary='pmemd.cuda'
    else:
        print 'Error: MD mode ' + MD_mode + ' not known.'

       
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
    
    os.system('echo ' + exec_amber_string + ' >> ' + work_dir + 'debug/MD_AMBER_command_line.log')
    os.system(exec_amber_string)
    
    if debug==False:
        os.remove(dir_trajectory_coords)
        os.remove(dir_outfile)    
        
def calc_coord(work_dir,MD_mode,debug,trajectory_index,trajectory_parent_index):
    pass       