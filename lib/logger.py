import pickle
import shutil
import os 
import re 
import glob
import lib.iteration as iteration

class Logger():
    # Global class variables for log filename formatting
    ITERATION_FNAME_FORMAT   = "{dir}/{id:05d}.iter"
    ITERATION_FNAME_WILDCARD = "{dir}/*.iter"
    ITERATION_FNAME_REGEXP   = ".*?(\d{5})\.iter"
    CONFIG_FNAME_FORMAT      = "{dir}/{id:05d}.conf"
    CONFIG_FNAME_REGEXP      = ".*?(\d{5})\.conf"
    
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
        # Sanitize the iteration_id if necessary
        if iteration_id < 0:
            iteration_id = self.getLastIterationId()
        iteration_fname = self.ITERATION_FNAME_FORMAT.format(                                 
                               dir=self.logdir,
                               id=iteration_id)
        return self.__loadIterationFile(iteration_fname)

    def loadConfigFile(self, iteration_id):
        """
        @return Filename of a config file for iteration_id
        """
        # Sanitize the iteration_id if necessary
        if iteration_id < 0:
            iteration_id = self.getLastIterationId()
        return self.CONFIG_FNAME_FORMAT.format(dir=self.logdir, id=iteration_id)
        
    def loadLastIterations(self, N=1):
        """
        @return List of last N iterations
        """
        last_iteration_id = self.getLastIterationId()
        return self.loadIterations(last_iteration_id - N + 1, last_iteration_id)
    
    def loadIterations(self, begin=0, end=-1):
        """
        @return list of iterations
        """
        iterations = []
        iteration_id = begin
        #for only reading the last iteration without knowing the id of the last iteration
        if (begin==-1 and end==-1):
            iteration = self.loadIteration(iteration_id)
            iterations.append(iteration)
        else:            
            
            while (os.path.isfile(
                           self.ITERATION_FNAME_FORMAT.format(
                                    dir=self.logdir, 
                                    id=iteration_id))
                  and (iteration_id <= end or end==-1) ):
                iteration = self.loadIteration(iteration_id)
                iterations.append(iteration)
                iteration_id += 1
            
        # We want to raise an error if... 
        # no iterations were read in
        if len(iterations) == 0:
            raise Exception("Could not read any iteration files "\
                            "for desired range {b:05d}-{e:05d}!".format(b = begin,
                                                                        e = end))          
        # or the range is known but not all files present 
        if end > -1: 
            if len(iterations) != end - begin + 1:
                raise Exception("Could not read all iteration files "\
                                "for desired range {b:05d}-{e:05d}!\n"\
                                "           Found iterations {first_found:05d}-{last_found:05d} "\
                                "of desired range.".format(b = begin,
                                                           e = end,
                                                           first_found = iterations[0].getId(),
                                                           last_found = iterations[-1].getId()))
        
        return iterations

    def getLastIterationId(self):
        """
        Returns the id of the last iteration found in log folder
        """
        last_iteration_id = 0
        for iteration_fname in glob.glob(self.ITERATION_FNAME_WILDCARD.format(dir=self.logdir)):
            this_id = re.search(self.ITERATION_FNAME_REGEXP, iteration_fname).group(1)
            if this_id > last_iteration_id:
                last_iteration_id = this_id
        return int(last_iteration_id)


