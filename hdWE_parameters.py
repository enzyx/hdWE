import json

#############################
# parsing the commandline
   
class HdWEParameters():
    """
    container for hdWE runtime paramters
    """
    def __init(self):
        pass
    
    def loadParams(self, workdir,
                         md_package,
                         max_iterations,
                         segments_per_bin,
                         minimal_probability,
                         coordinate_threshold,
                         max_bins,
                         starting_structure,
                         convergence_check_range,
                         convergence_check_threshold,
                         logfile = "hdWE.log",
                         debug = False,
                         reweighting_range = 0):
        self.workdir              = str(workdir)                # str, DIR
        self.logfile              = str(logfile)                 # str, FILE
        self.max_iterations       = int(max_iterations)          # int
        self.coordinate_threshold = float(coordinate_threshold)  # float
        self.segments_per_bin     = int(segments_per_bin)        # int
        self.minimal_probability  = float(minimal_probability)   # float
        self.max_bins             = int(max_bins)                # int
        self.debug                = bool(debug)                  # bool
        self.md_package           = str(md_package)              # str
        self.reweighting_range    = int(reweighting_range)       # int
        self.convergence_check_range     = int(convergence_check_range)       # int
        self.convergence_check_threshold = int(convergence_check_threshold)   #float
        self.starting_structure   = str(starting_structure)      # str
        
    def loadConfParameters(self, config, debug=False):
        """
        load hdWE parameters from config file
        """
        # define md_package string (amber/gromacs)
        if "amber" in config.sections():
            md_package = "amber"
        elif "gromacs" in config.sections():
            md_package = "gromacs"
        else:
            raise Exception("No MD package (amber/gromacs) section in configuration file")
        self.md_package           = md_package
        self.workdir              = str(config.get('hdWE','workdir'))
        self.guaranteeWorkdirSlash()
        self.max_iterations       = int(config.get('hdWE','max-iterations'))
        self.segments_per_bin     = int(config.get('hdWE','segments-per-bin'))
        self.minimal_probability  = float(config.get('hdWE','minimal-probability'))
        self.coordinate_threshold = float(config.get('hdWE','threshold'))
        self.max_bins             = int(config.get('hdWE','max-bins'))
        self.logfile              = self.workdir + str(config.get('hdWE','logfile'))
        self.debug                = debug
        self.reweighting_range    = int(config.get('hdWE','reweighting-range'))
        self.convergence_check_range     = int(config.get('hdWE','convergence-check-range'))       
        self.convergence_check_threshold = float(config.get('hdWE','convergence-check-threshold'))   
        self.starting_structure   = str(config.get('hdWE','starting-structure'))
        
    def loadJsonParams(self, json_string):
        param_dict = json.loads(json_string)
        self.workdir              = param_dict.get("workdir")                # str, DIR
        self.guaranteeWorkdirSlash()
        self.logfile              = self.workdir + param_dict.get("logfile")  # str, FILE
        self.max_iterations       = param_dict.get("max_iterations")          # int
        self.coordinate_threshold = param_dict.get("coordinate_threshold")    # float
        self.segments_per_bin     = param_dict.get("segments_per_bin")        # int
        self.minimal_probability  = param_dict.get("minimal_probability")     # float
        self.max_bins             = param_dict.get("max_bins")                # int
        self.debug                = param_dict.get("debug")                   # bool
        self.md_package           = param_dict.get("md_package")              # str
        self.reweighting_range    = param_dict.get("reweighting_range")       # int
        self.starting_structure   = param_dict.get("starting_structure")      # str
        self.convergence_check_range     = param_dict.get("convergence-check-range")     # int      
        self.convergence_check_threshold = param_dict.get("convergence-check-threshold") # float
        
    def getLogString(self):
        return json.dumps(self.__dict__, sort_keys=True)
        
    def guaranteeWorkdirSlash(self):
        """
        guarantee a working workdir variable
        """
        if self.workdir and self.workdir[-1] != "/":
            self.workdir +="/"
