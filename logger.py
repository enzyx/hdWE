#!/usr/bin/python3
from __future__ import print_function
import sys, os
import glob
import json
from segment import Segment
from bin import Bin
from iteration import Iteration
from hdWE_parameters import HdWEParameters

class Logger():
    """
        writes iterations to file and reads them.
        
        HOWTO:
        # start a logger:
        logger = Logger(logfilename = "logfile.log")
        
        # log a set of command line arguments
        logger.logParameters(args)
        
        # log an array of iterations
        logger.logIterations(iterations)
        
        # log a single iteration
        logger.logIteration(iteration)
        
        # load parameters from the logfile:
        logger.loadHdWEParameters(args)
        
        # load (correctly) logged iterations
        read_iterations = logger.load_iterations()
        
        # close the logger(!):
        logger.close()
    
    """
    def __init__(self, filename, writeable_logfile = False, append = True, debug = False):
        self.logfilename = filename
        self.debug = debug

        # open logfile
        if create_logfile:
            self.logfile = open(self.logfilename, "a+")
        else:
            self.logfile = open(self.logfilename, "r")
    
    def logParameters(self, hdWE_parameters):
        paramline = hdWE_parameters.getLogString()
        self.logfile.write("hdWE_parameters|"+paramline + "\n")
        self.logfile.flush()
            
    def logIteration(self, iteration):
        """
        appends a single iteration to the logfile
        """
        self.logfile.write("Iteration {it:05d} - {nbins} bins\n".format(\
            it    = iteration.getId(),
            nbins = iteration.getNumberOfBins()))
        for _bin in iteration.bins:
            self._logBin(_bin)
        self.logfile.write("completed iteration {it:05d}\n".format(\
                            it = iteration.getId()))
        
    def logIterations(self, iterations):
        """
        writes all interations to the logfile
        """
        for iteration in iterations:
            self.logIteration(iteration)
    
    def loadHdWEParameters(self):
        """
        @return hdWEParameters instance with parameters from logfile
        """
        hdWE_parameters = HdWEParameters()
        hdWE_parameters.loadLogfileParameters(self._getHdWEParameterString())
        return hdWE_parameters
            
    def loadIterations(self, first=0, last=-1, bCheckFiles=True):
        """
            @return read-in iterations
        """
        sys.stdout.write("reading iterations from {file}\n".format(file = self.logfilename))
        iterations = []
        iteration = Iteration(-1)
        number_of_iterations_read = 0   # number of iterations read
        b_read_iteration = False        # True if iteration is to be read
        
        # read in hdWE parameters (e.g. for md module for file testing)
        self._createOwnHdWEParameters()
        
        with open(self.logfilename, "r") as readfile:
            file_lines = readfile.readlines()
            # compare arguments, suspended for now because the logger
            # normally doesn't know all arguments
            #~ argsline = file_lines[0]
            #~ self.__checkArguments(json.loads(argsline))

            for line in file_lines:
                #~ print ("len(iterations):",len(iterations))
                #~ print("number of iterations read: ", number_of_iterations_read)
                if line.strip().startswith(("@","#")):
                    continue
                if line.strip()== "":
                    continue
                if "hdWE_parameters" in line:
                    continue
                
                # only load itererations if completed flag exists
                if line[0:9].lower() == "completed" and line.split()[-1] == "{id:05d}".format(id=iteration.getId()):
                    if self.__checkIteration(iteration,
                                             target_number_of_read_bins,
                                             first,
                                             last,
                                             bCheckFiles):
                        iterations.append(iteration)
                        
                # read iteration
                if line[0:9].lower() == "iteration":
                    # parse iteration line    
                    iteration_line = line.split()
                    iteration=Iteration(int(iteration_line[1]))
                    target_number_of_read_bins = int(iteration_line[3])
                    number_of_iterations_read += 1
                    # check if iteration is in range
                    if iteration.getId() >= first and \
                       (last == -1 or iteration.getId() <= last):
                        b_read_iteration = True
                    else:
                        b_read_iteration = False
                # read bin data
                if b_read_iteration and line[:5] == "{\"b\":":
                    read_bin = self._loadBin(line)
                    # check for iteration_id consistency
                    if self.debug and read_bin.getIterationId() != iteration.getId():
                        raise Exception("Iteration mismatch: Iteration: {it_id},"+
                                        "Bin: {bin_it_it}".format(\
                                            it_id=self.iteration.getId(),
                                            bin_it_id = read_bin.getIterationId()))
                    # append bin
                    iteration.bins.append(read_bin)

        if self.debug:
            sys.stdout.write("finished reading iterations\n")    
        return iterations
     
    def _getHdWEParameterString(self):
        """
        @return read in arguments line
        """
        # read in arguments from logfile
        for line in reversed(open(self.logfilename).readlines()):
            if "hdWE_parameters" in line:
                param_line = line.split("|")[1]
                break
                
        return param_line
    
    def _createOwnHdWEParameters(self):
        self.hdWE_parameters = self.loadHdWEParameters()
        if self.hdWE_parameters.md_package.lower() == "amber":
            from amber_module import MD_analysis_module
            self.md_analysis_module = MD_analysis_module(
                                        workdir = self.hdWE_parameters.workdir,
                                        jobname = self.hdWE_parameters.jobname,
                                        debug = self.debug)
        if self.hdWE_parameters.md_package.lower() == "gromacs": 
            pass       
       
    def _logBin(self, _bin):
        """
            writes a bin line to the logfile
        """
        self.bin_part = json.dumps(_bin, default=self._convertBin,
                             sort_keys=True, separators=(',',':'))
        self.segments_part = ""
        for segment in _bin:
            self.segments_part = json.dumps(_bin.segments,
                                   default=self._convertSegment,
                                   sort_keys=True, separators=(',',':'))

        line = self.bin_part+"|"+self.segments_part+"\n"
        self.logfile.write(line)
        self.logfile.flush()     
    
    def _loadBin(self, line):
        """
            reads a given bin line
        """
        splitline = line.split("|")
        
        # read bin data
        bin_string = splitline[0]
        self.newbin = json.loads(bin_string, object_hook=self._reconvertBin)
        
        
        # read segment data
        if splitline[1][0] == "[":
            segments_string = splitline[1]
            segment_list = json.loads(segments_string)
            for segment_dictionary in segment_list:        
                if set(('p', 'i', 'b', 's')) <= set(segment_dictionary):
    
                    # check for consistency: note that segment_dictionary.get('i') is the parent iteration!
                    #~ if int(self.newbin.getIterationId()) != int(segment_dictionary.get('i')+1):
                        #~ raise Exception("Iteration_Ids of segment ({seg}) and bin ({bin}) mismatch!".\
                                        #~ format(seg=int(segment_dictionary.get('i')+1),
                                        #~ bin=int(self.newbin.getIterationId())))
    
                    self.newbin.generateSegment(\
                        probability = segment_dictionary.get('p'),
                        parent_bin_id = segment_dictionary.get('b'),
                        parent_segment_id = segment_dictionary.get('s'))            
                else:
                    raise Exception("Non-segment or insufficient data parsed as segment:\n" + segment_dictionary)
                
        return self.newbin
        
    def _convertSegment(self, segment):
        """
            converts parent data and probability to built-in format
            @return dictionary of segment data to save
        """
        self.dictionary = {'p':segment.getProbability(),
            'i':segment.getParentIterationId(),
            'b':segment.getParentBinId(),
            's':segment.getParentSegmentId()}
        return self.dictionary
    
    def _reconvertSegment(self, dictionary, iteration_id, bin_id, segment_id):
        """
            reads stored segment data
            @return segment from saved data
        """
        if set(('p', 'i', 'b', 's')) <= set(dictionary):
            if iteration_id != dictionary.get('i')+1:
                raise Exception("Iteartion ID mismatch: Segment: {seg} != Bin: {_bin}".\
                                 format(seg=dictionary.get('i')+1, _bin=iteration_id))
            segment = Segment(probability=float(dictionary.get('p')),
                              parent_bin_id=int(dictionary.get('b')),
                              parent_segment_id=int(dictionary.get('s')),
                              iteration_id=iteration_id,
                              bin_id=bin_id,
                              segment_id=segment_id)
        else:
            raise Exception("{module}: No Segment given for reconversion.".\
                             format(module=self))
        return segment

    def _convertBin(self, _bin):
        """
        Convert bin data (except segments) to built-in types
        @return list of dictionaries with bindata
        """
        self.bin_dictionary = {\
            'i':_bin.getIterationId(),
            'b':_bin.getId(),
            'r':_bin.getReferenceNameString(),
            'n':_bin.getTargetNumberOfSegments(),
            'conv':_bin.isConverged()}
        
        return self.bin_dictionary
    
    def _reconvertBin(self, bin_dictionary):
        """
            @return returns a bin without segments from saved data
        """
        # read bin data
        if set(('i', 'b', 'r', 'n')) <= set(bin_dictionary):
            self.reference_name           = bin_dictionary.get('r').split("_")
            self.newbin = Bin(\
                iteration_id              = int(bin_dictionary.get('i')), 
                bin_id                    = int(bin_dictionary.get('b')), 
                reference_iteration_id    = int(self.reference_name[0]), 
                reference_bin_id          = int(self.reference_name[1]),
                reference_segment_id      = int(self.reference_name[2]),
                target_number_of_segments = int(bin_dictionary.get('n')),
                outrates_converged        = bool(bin_dictionary.get('conv')) )
        else:
            raise Exception("Non-bin or insufficient data parsed as bin.")

        return self.newbin
        
    def printIteration(self, iteration):
        """
            prints a human (barely) readable iteration
        """
        for _bin in iteration.bins:
            self.printBin(_bin)     

    def printBin(self, _bin):
        """
            prints bin data
        """
        print ("____Bin:____")
        print ("Bin {i:05d}_{b:05d}".format(i=_bin.getIterationId(), b=_bin.getId()))
        print ("Reference: {ref}".format(ref=_bin.getReferenceNameString()))
        print ("Target number of segments: {t}".format(t=_bin.getTargetNumberOfSegments()))
        print ("____Segments____:")
        for segment in _bin.segments:
            self.print_segment(segment)            

    def print_segment(self, segment):
        """
            prints segment data
        """
        print ( "segment:", segment.getNameString() )
        print ( "probability =", segment.getProbability() )
        print ( "parent =", segment.getParentNameString() )
                
    def close(self):
        """
            closes the logfile
        """
        self.logfile.close()
        
    def __checkArguments(self, read_args):
        """
        checks if the arguments of the logfile match the input arguments
        @ return True if arguments match
        """
        bCheck = True
        if self.args.amber != read_args.get("amber") or self.args.gromacs != read_args.get("gromacs"):
            # find out what the arguments and the logfile want to use
            if self.args.amber:
                arg_package = "amber"
            elif self.args.gromacs:
                arg_package = "gromacs"
            else:
                arg_package = "none"
            if read_args.get("amber") == "true":
                log_package = "amber"
            elif read_args.get("gromacs") == "true":
                log_package = "gromacs"
            else:
                log_package = "none"            
            sys.stderr.write("WARNING: Simulation package (amber/gromacs)" +
                                " of logfile ({log})".format(\
                                log = log_package) +
                                " and arguments ({arg}) don't match.\n".format(\
                                arg = arg_package))
            bCheck = False                  
        if self.args.coordinate_threshold != read_args.get("coordinate_threshold"):
            sys.stderr.write("WARNING: threshold" +
                                " of logfile ({log}) and arguments ({arg}) don't match.\n".format(\
                                log = read_args.get("coordinate_threshold"),
                                arg = self.args.coordinate_threshold ))                
            bCheck = False                
        if self.args.input_minimal_probability != read_args.get("input_minimal_probability"):
            sys.stderr.write("WARNING: minimal probability" +
                                " of logfile ({log}) and arguments ({arg}) don't match.\n".format(\
                                log = read_args.get("input_minimal_probability"),
                                arg =  self.args.input_minimal_probability ))                
            bCheck = False                
        if self.args.segments_per_bin != read_args.get("segments_per_bin"):
            sys.stderr.write("WARNING: segments_per_bin" +
                                " of logfile ({log}) and arguments ({arg}) don't match.\n".format(\
                                log = read_args.get("segments_per_bin"),
                                arg =  self.args.segments_per_bin ))                
            #~ bCheck = False
            
        return bCheck
    def __checkIteration(self, iteration, target_number_of_read_bins, first, last, bCheckFiles):
        """
        @return True if iteration is complete (all bins and 1.0 Probability)
        """
        bNbins = False
        bProbability = True
        bRange = True
        # check if iteration is in range
        if iteration.getId() < first or (last != -1 and iteration.getId() > last):
                bRange = False
                return False
        # check iteration number of bins
        if iteration.getNumberOfBins() == target_number_of_read_bins:
            bNbins = True
        else:
            sys.stderr.write("read-in error: Bin number mismatch: Iteration {it} expected: {target},"\
                             " read {nbins} bins\n".format(
                                target = target_number_of_read_bins,
                                it = iteration.getId(),
                                nbins = iteration.getNumberOfBins()))
        # check iteration probability
        if self.debug:
            if iteration.checkProbability():
                bProbability = True
            else:
                sys.stderr.write("read-in error: Wrong total Iteration probability: {prob}\n".format(\
                                prob = iteration.getProbability())) 
        # check iteration segment files
        bAllFilesPresent = True
        if bCheckFiles:
            for segment in iteration.getSegments():
                    if not self.md_analysis_module.isSegmentFile(segment):
                        bAllFilesPresent = False
                        sys.stderr.write("read-in error: File missing for segment {}\n".format(\
                                    segment.getNameString()))
                        break
        #~ # append iteration
        if bNbins and bProbability and bRange and bAllFilesPresent:
            return True
        else:
            sys.stderr.write("omitting iteration {it}\n".format(it=iteration.getId()))
        return False
