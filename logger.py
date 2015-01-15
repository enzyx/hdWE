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
        self.iterations = []
        with open(self.logfilename, "r") as readfile:
            # find first iteration number from first bin:
            self.firstline = readfile.readline()
            self.read_bin = self._load_bin(self.firstline)
            self.iteration_id = self.read_bin.getIterationId()
            self.iteration = Iteration(self.iteration_id)
            readfile.seek(0)
            for self.line in readfile:
                """
                 for now json basically just starts self._reconvert_bin where a
                 self.newbin is created
                """
                self.read_bin = self._load_bin(self.line)
                                
                # create new iteration if bin_id changed
                if self.read_bin.getIterationId() != self.iteration_id:
                    self.iterations.append(self.iteration)
                    self.iteration_id = self.read_bin.getIterationId()
                    self.iteration = Iteration(self.iteration_id)
                    
                self.iteration.bins.append(self.read_bin)
            self.iterations.append(self.iteration)
        return self.iterations
       
    def _log_bin(self, _bin):
        self.bin_part = json.dumps(_bin, default=self._convert_bin,
                             sort_keys=True, separators=(',',':'))
        self.segments_part = ""
        for segment in _bin:
            self.segments_part = json.dumps(_bin.segments,
                                   default=self._convert_segment,
                                   sort_keys=True, separators=(',',':'))

        self.line = self.bin_part+"|"+self.segments_part+"\n"
        self.logfile.write(self.line)
        self.logfile.flush()     
    
    def _load_bin(self, line):
        sline = line.split("|")
        bin_string = sline[0]
        segments_string = sline[1]
        
        # read bin data
        self.newbin = json.loads(bin_string, object_hook=self._reconvert_bin)
        
        
        # read segment data
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
            @return dictionary of data to save
        """
        self.dictionary = {'p':segment.getProbability(),
            'i':segment.getParentIterationId(),
            'b':segment.getParentBinId(),
            's':segment.getParentSegmentId()}
        return self.dictionary
    
    def _reconvert_segment(self, dictionary, iteration_id, bin_id, segment_id):
        """
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
        Convert bin to built-in types
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
            @return returns a bin from saved data
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
        for _bin in iteration.bins:
            self.print_bin(_bin)     

    def print_bin(self, _bin):
        print ("____Bin:____")
        print ("Bin {i:05d}_{b:05d}".format(i=_bin.getIterationId(), b=_bin.getId()))
        print ("Reference: {ref}".format(ref=_bin.getReferenceNameString()))
        print ("Target number of segments: {t}".format(t=_bin.getTargetNumberOfSegments()))
        print ("____Segments____:")
        for segment in _bin.segments:
            self.print_segment(segment)            

    def print_segment(self, segment):
        print ( "segment:", segment.getNameString() )
        print ( "probability =", segment.getProbability() )
        print ( "parent =", segment.getParentNameString() )
                
    def close(self):
        self.logfile.close()
