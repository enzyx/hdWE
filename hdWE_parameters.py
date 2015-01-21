import json

#############################
# parsing the commandline

class HdWEParameters():
    """
    container for hdWE runtime paramters
    """
    def __init(self):
        pass
    
    def loadParams(self, work_dir,
                        md_conf,
                        md_package,
                        max_iterations,
                        segments_per_bin,
                        minimal_probability,
                        coordinate_threshold,
                        max_bins,
                        logfile = "hdWE.log",
                        debug = False):
        self.work_dir             = str(work_dir)                # str, DIR
        self.md_conf              = str(md_conf)                 # str, FILE
        self.logfile              = str(logfile)                 # str, FILE
        self.max_iterations       = int(max_iterations)          # int
        self.coordinate_threshold = float(coordinate_threshold)  # float
        self.segments_per_bin     = int(segments_per_bin)        # int
        self.minimal_probability  = float(minimal_probability)   # float
        self.max_bins             = int(max_bins)                # int
        self.debug                = bool(debug)                  # bool
        self.md_package           = str(md_package)              # str
        
    def loadJsonParams(self, json_string):
        param_dict = json.loads(json_string)
        self.work_dir             = param_dict.get("work_dir")                # str, DIR
        self.md_conf              = param_dict.get("md_conf")                 # str, FILE
        self.logfile              = param_dict.get("logfile")                 # str, FILE
        self.max_iterations       = param_dict.get("max_iterations")          # int
        self.coordinate_threshold = param_dict.get("coordinate_threshold")    # float
        self.segments_per_bin     = param_dict.get("segments_per_bin")        # int
        self.minimal_probability  = param_dict.get("minimal_probability")     # float
        self.max_bins             = param_dict.get("max_bins")                # int
        self.debug                = param_dict.get("debug")                   # bool
        self.md_package           = param_dict.get("md_package")              # str

    def getLogString(self):
        return json.dumps(self.__dict__, sort_keys=True)
