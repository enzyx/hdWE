class propagate:
    """trajectory propagation via MD."""
    
    
    def __init__(self,work_dir,MD_software='amber',MD_mode='pmemd',MD_n_parallel=1,MD_debug=False):
        """sets general propagation variables"""  
        
        self.work_dir       =   work_dir         # working directory
        self.MD_mode        =   MD_mode          # specifies whether to run the MD with cpu or cuda or ...
        self.MD_n_parallel  =   MD_n_parallel       # speficies how many MD simulations to run in parallel 
        self.MD_software    =   MD_software      # specifies which MD software is used.
        self.MD_debug       =   MD_debug            # debug flag

    def run_MD(self,trajectory):
        """runs the MD simulations."""

        #Load the MD module         
        if self.MD_software=='amber':
            import amber_module as MD_module
        else:
            print 'Error: MD software ' + self.MD_software + ' not known.' 

                    
        if self.n_parallel==1:
	    for Bin_index in range(0,len(Bin)):
            	for Tra_index in range(0,len(Bin[Bin_index].Trajectory)):
                	MD_module.runMD(self.work_dir,self.MD_mode,self.debug, \ 
					iteration,Bin_index,Tra_index,Bin[Bin_index].Trajectory[Tra_index].parent_bin,Bin[Bin_index].Trajectory[Tra_index].parent_trajectory)
                	#here: print status of propagation to log file

        elif  self.n_parallel>1:     
            print 'Not yet implemented'                

        else: 
            print 'Error with n_parallel'                













    

    
                            
    
    
    
