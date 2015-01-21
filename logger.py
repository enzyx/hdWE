#!/usr/bin/python3

import sys, os
import json
from segment import Segment
from bin import Bin
from iteration import Iteration


class Logger():
    """
        writes iterations to file and reads them.
        
        HOWTO:
        # start a logger:
        logger = Logger(logfilename = "logfile.log")
        
        # log a set of command line arguments
        logger.log_arguments(args)
        
        # log an array of iterations
        logger.log_iterations(iterations)
        
        # log a single iteration
        logger.log_iteration(iteration)
        
        # load arguments from the logfile:
        logger.load_arguments(args)
        
        # load (correctly) logged iterations
        read_iterations = logger.load_iterations()
        
        # close the logger(!):
        logger.close()
    
    """
    def __init__(self, filename = "logfile.log"):
        self.logfilename = filename
        self.logfile = open(self.logfilename, "a+")
    
    def log_arguments(self, args):
        argline = json.dumps(vars(args), self.logfile, sort_keys=True)
        self.logfile.write(argline + "\n")
        self.logfile.flush()
            
    def log_iteration(self, iteration):
        """
        appends a single iteration to the logfile
        """
        self.logfile.write("Iteration {it:05d} - {nbins} bins\n".format(\
            it    = iteration.getId(),
            nbins = iteration.getNumberOfBins()))
        for _bin in iteration.bins:
            self._log_bin(_bin)
        
    def log_iterations(self, iterations):
        """
        writes all interations within [first, last] to the logfile
        """
        for iteration in iterations:
            self.log_iteration(iteration)
            
    def load_arguments(self, args):
        """
        @return read in arguments
        """
        # read in crucial arguments from logfile
        for line in reversed(open(self.logfilename).readlines()):
            if "coordinate_threshold" in line:
                read_args = json.loads(line)
                break
        
        # tell the user what we're doing
        sys.stderr.write("WARNING: overwriting input arguments with last logged ones:\n")
        sys.stderr.write("amber = {}\n".format(str(read_args.get("amber"))))
        sys.stderr.write("gromacs = {}\n".format(read_args.get("gromacs")))
        sys.stderr.write("coordinate_threshold = {}\n".format(read_args.get("coordinate_threshold")))
        sys.stderr.write("input_minimal_probability = {}\n".format(read_args.get("input_minimal_probability")))
        sys.stderr.write("segments_per_bin = {}\n".format(read_args.get("segments_per_bin")))
        sys.stderr.write("input_md_conf = {}\n".format(read_args.get("input_md_conf")))
        sys.stderr.write("work_dir = {}\n".format(read_args.get("work_dir")))
        
        args.amber = read_args.get("amber")
        args.gromacs = read_args.get("gromacs")
        args.coordinate_threshold = read_args.get("coordinate_threshold")
        args.input_minimal_probability = read_args.get("input_minimal_probability")
        args.segments_per_bin = read_args.get("segments_per_bin") 
        args.input_md_conf = read_args.get("input_md_conf")
        args.work_dir = read_args.get("work_dir")
                
        return args
            
    def load_iterations(self, first=0, last=-1):
        """
            @return read-in iterations
        """
        iterations = []
        iteration = Iteration(-1)
        number_of_iterations_read = 0
        with open(self.logfilename, "r") as readfile:
            file_lines = readfile.readlines()
            # compare arguments, suspended for now because the logger
            # normally doesn't know all arguments
            #~ argsline = file_lines[0]
            #~ self.__check_arguments(json.loads(argsline))

            for line in file_lines:
                #~ print ("len(iterations):",len(iterations))
                #~ print("number of iterations read: ", number_of_iterations_read)
                if line.strip().startswith(("@","#")):
                    continue
                if line.strip()== "":
                    continue
                if "coordinate_threshold" in line:
                    continue
                
                # read iterations
                if line[0:9].lower() == "iteration":
                    # check and append previous iteration if it's not the first one
                    if number_of_iterations_read > 0:
                        if self.__check_iteration(iteration,
                                                  target_number_of_read_bins,
                                                  first,
                                                  last):
                            iterations.append(iteration)
                    # parse iteration line    
                    iteration_line = line.split()
                    iteration=Iteration(int(iteration_line[1]))
                    target_number_of_read_bins = int(iteration_line[3])
                    number_of_iterations_read += 1

                # read bin data
                if line[:5] == "{\"b\":":
                    read_bin = self._load_bin(line)
                    # check for iteration_id consistency
                    if read_bin.getIterationId() != iteration.getId():
                        raise Exception("Iteration mismatch: Iteration: {it_id},"+
                                        "Bin: {bin_it_it}".format(\
                                            it_id=self.iteration.getId(),
                                            bin_it_id = read_bin.getIterationId()))
                    # append bin
                    iteration.bins.append(read_bin)

            # append last iteration
            if self.__check_iteration(iteration, target_number_of_read_bins, first, last):
                iterations.append(iteration)
                
        return iterations
       
    def _log_bin(self, _bin):
        """
            writes a bin line to the logfile
        """
        self.bin_part = json.dumps(_bin, default=self._convert_bin,
                             sort_keys=True, separators=(',',':'))
        self.segments_part = ""
        for segment in _bin:
            self.segments_part = json.dumps(_bin.segments,
                                   default=self._convert_segment,
                                   sort_keys=True, separators=(',',':'))

        line = self.bin_part+"|"+self.segments_part+"\n"
        self.logfile.write(line)
        self.logfile.flush()     
    
    def _load_bin(self, line):
        """
            reads a given bin line
        """
        splitline = line.split("|")
        
        # read bin data
        bin_string = splitline[0]
        self.newbin = json.loads(bin_string, object_hook=self._reconvert_bin)
        
        
        # read segment data
        if splitline[1][0] == "[":
            segments_string = splitline[1]
            segment_list = json.loads(segments_string)
            for segment_dictionary in segment_list:        
                if set(('p', 'i', 'b', 's')) <= set(segment_dictionary):
    
                    # check for consistency: note that segment_dictionary.get('i') is the parent iteration!
                    if int(self.newbin.getIterationId()) != int(segment_dictionary.get('i')+1):
                        raise Exception("Iteration_Ids of segment ({seg}) and bin ({bin}) mismatch!".\
                                        format(seg=int(segment_dictionary.get('i')+1),
                                        bin=int(self.newbin.getIterationId())))
    
                    self.newbin.generateSegment(\
                        probability = segment_dictionary.get('p'),
                        parent_bin_id = segment_dictionary.get('b'),
                        parent_segment_id = segment_dictionary.get('s'))            
                else:
                    raise Exception("Non-segment or insufficient data parsed as segment:\n" + segment_dictionary)
                
        return self.newbin
        
    def _convert_segment(self, segment):
        """
            converts parent data and probability to built-in format
            @return dictionary of segment data to save
        """
        self.dictionary = {'p':segment.getProbability(),
            'i':segment.getParentIterationId(),
            'b':segment.getParentBinId(),
            's':segment.getParentSegmentId()}
        return self.dictionary
    
    def _reconvert_segment(self, dictionary, iteration_id, bin_id, segment_id):
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

    def _convert_bin(self, _bin):
        """
        Convert bin data (except segments) to built-in types
        @return list of dictionaries with bindata
        """
        self.bin_dictionary = {\
            'i':_bin.getIterationId(),
            'b':_bin.getId(),
            'r':_bin.getReferenceNameString(),
            'n':_bin.getTargetNumberOfSegments()}
        
        return self.bin_dictionary
    
    def _reconvert_bin(self, bin_dictionary):
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
                target_number_of_segments = int(bin_dictionary.get('n')) )
        else:
            raise Exception("Non-bin or insufficient data parsed as bin.")

        return self.newbin
        
    def print_iteration(self, iteration):
        """
            prints a human (barely) readable iteration
        """
        for _bin in iteration.bins:
            self.print_bin(_bin)     

    def print_bin(self, _bin):
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
        
    def __check_arguments(self, read_args):
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
    def __check_iteration(self, iteration, target_number_of_read_bins, first, last):
        """
        @return True if iteration is complete (all bins and 1.0 Probability)
        """
        bNbins = False
        bProbability = True
        bRange = True
        # check iteration number of bins
        if iteration.getNumberOfBins() == target_number_of_read_bins:
            bNbins = True
        else:
            raise Exception("Bin number mismatch: Iteration had: {target},".format(\
                        target = target_number_of_read_bins)+
                        " read {nbins} bins".format(\
                            nbins = iteration.getNumberOfBins()))
        # check iteration probability
        if iteration.checkProbability():
            bProbability = True
        else:
            sys.stderr.write("Wrong total Iteration probability: {prob}\n".format(\
                            prob = iteration.getProbability())) 
        # check if iteration is in range
        if iteration.getId() < first:
                bRange = False
        elif last != -1 and iteration.getId() > last:
                bRange = False
        # append iteration
        if bNbins and bProbability and bRange:
            return True
        return False
