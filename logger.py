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
        logger = Logger(logfilename = "test.log", read_filename = "test.log")
        
        # log an array of iterations
        logger.log_iterations(iterations)
        
        # log a single iteration
        logger.log_iteration(iteration)
        
        # read a (correctly) logged file
        read_iterations = logger.load_iterations()
        
        # close the logger(!):
        logger.close()
    
    """
    def __init__(self, filename):
        self.logfilename = filename
        self.logfile = open(self.logfilename, "a+")
            
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
        writes all interations to the logfile
        """
        for iteration in iterations:
            self.log_iteration(iteration)
            
    def load_iterations(self):
        """
            @return read-in iteration
        """
        def check_iteration(iteration, target_number_of_read_bins):
            """
            @return True if iteration is complete (all bins and 1.0 Probability)
            """
            bNbins = False
            bProbability = False
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
                raise Exception("Wrong total Iteration probability: {prob}".format(\
                                prob = iteration.getProbability())) 
            # append iteration
            if bNbins and bProbability:
                return True
            return False
        
        
        self.iterations = []
        with open(self.logfilename, "r") as readfile:
            for line in readfile:
				if line.strip() == "":
					continue
                if line[0:9].lower() == "iteration":
                    # check and append previous iteration if it exists
                    try:
                        self.iteration
                    except:
                        pass
                    else:
                        if check_iteration(self.iteration, self.target_number_of_read_bins):
                            self.iterations.append(self.iteration)
                    # parse iteration line    
                    self.iteration_line = line.split()
                    self.iteration=Iteration(int(self.iteration_line[1]))
                    self.target_number_of_read_bins = int(self.iteration_line[3])

                else:
                    # read bin data
                    self.read_bin = self._load_bin(line)
                    # check for iteration_id consistency
                    if self.read_bin.getIterationId() != self.iteration.getId():
                        raise Exception("Iteration mismatch: Iteration: {it_id},"+
                                        "Bin: {bin_it_it}".format(\
                                            it_id=self.iteration.getId(),
                                            bin_it_id = self.read_bin.getIterationId()))
                    # append bin
                    self.iteration.bins.append(self.read_bin)

            # append last iteration
            if check_iteration(self.iteration, self.target_number_of_read_bins):
                self.iterations.append(self.iteration)
            
        return self.iterations
       
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
