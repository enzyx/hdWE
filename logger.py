import pickle
import shutil
import os 
import re 
import glob

class Logger():
    # Global class variables for log filename formatting
    ITERATION_FNAME_FORMAT = "{dir}/{id:05d}.iter"
    ITERATION_FNAME_REGEXP = ".*?(\d{5})\.iter"
    CONFIG_FNAME_FORMAT    = "{dir}/{id:05d}.conf"
    CONFIG_FNAME_REGEXP    = ".*?(\d{5})\.conf"
    
    def __init__(self, LOGDIR):
        """
        
        """
        self.logdir = "{ld}/".format(ld=LOGDIR)
        
    def logConfigFile(self, CONFIGFILE, iteration_id):
        """
        
        """
        shutil.copyfile(CONFIGFILE, self.CONFIG_FNAME_FORMAT.format(
                                 dir=self.logdir,
                                 id=iteration_id))
    
    def logIteration(self, iteration):
        """
        Save the iteration to an pickle object file
        """
        iteration_file = open(self.ITERATION_FNAME_FORMAT.format(
                                 dir=self.logdir,
                                 id=iteration.getId()), 
                              'w')
        pickle.dump(iteration, iteration_file, protocol=2)
        iteration_file.close()
    
    def log(self, iteration, CONFIGFILE):
        """
        Logs all necessary files:
            - Iteration
            - ConfigFile
        """
        self.logIteration(iteration)
        self.logConfigFile(CONFIGFILE, iteration.getId())
    
    def __loadIterationFile(self, fname):
        """
        Load the local copy of iteration from the 
        given file name 
        """
        iteration_file = open(fname, 'r')
        iteration = pickle.load(iteration_file)
        iteration_file.close()
        return iteration

    def loadIteration(self, iteration_id):
        """
        Load the local copy of iteration 
        for the given ID 
        """
        iteration_fname = self.ITERATION_FNAME_FORMAT.format(                                 
                               dir=self.logdir,
                               id=iteration_id)
        return self.__loadIterationFile(iteration_fname)

    def loadConfigFile(self, iteration_id):
        """
        """
        pass
    
    def loadLastIterations(self, N=1):
        """
        @return List of last N iterations
        """
        last_iteration_id = 0
        for iteration_fname in glob.glob(self.ITERATION_FNAME_FORMAT.format(dir=self.logdir,id="*")):
            this_id = re.search(self.ITERATION_FNAME_REGEXP, iteration_fname).group(1)
            if this_id > last_iteration_id:
                last_iteration_id = this_id
        return self.loadIterations(last_iteration_id - N + 1, last_iteration_id)
    
    def loadIterations(self, begin=0, end=-1):
        """
        @return list of iterations
        """
        iterations = []
        iteration_id = begin
        while (os.path.isfile(
                       self.ITERATION_FNAME_FORMAT.format(
                                dir=self.logdir, 
                                id=iteration_id))
              and (iteration_id <= end or end==-1) ):
            iteration = self.loadIteration(iteration_id)
            iterations.append(iteration)
            iteration_id += 1
        
        # We want to raise an error if 
        # the end is known but not all files present 
        if end < -1:
            if len(iterations) != end - begin:
                raise Exception("Could not read all iteration files "\
                                "for desired range {0:05d}-{1:05d}!\n"\
                                "Found first {0:d} iterations.".format(begin, end, len(iterations)))
        
        return iterations
