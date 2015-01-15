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
        logger = Logger(write_filename = "test.log", read_filename = "test.log")
        
        # log an array of iterations
        logger.log_iterations(iterations)
        
        # read a (correctly) logged file
        read_iterations = logger.read_iterations()
        
        # close the logger(!):
        logger.close()
    
    """
    def __init__(self, write_filename = "hdWE.log", read_filename = "hdWE.log"):
        self.write_filename = write_filename
        self.read_filename = read_filename
        self.logfile = open(self.write_filename, "a+", encoding='utf-8')
        #~ if os.path.isfile(self.read_filename):
            #~ self.readfile = open(self.read_filename, "r")

    
    def log_iteration(iteration):
        """
        writes a single iteration to the logfile
        """
        self.logfile = open(self.write_filename,"a")
        json.dump(iteration, self.logfile)
        self.logfile.flush()
        
    def log_all_iterations(self, iteration_array):
        """
        writes all interations to the logfile
        """
        for iteration in iteration_array:
            log_iteration(iteration)

    def load_iterations(self):
        """
            returns a read in iteration
        """
        iteration = Iteration()
        return iteration

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
            return segment
        else:
            raise Exception("{module}: No Segment given for reconversion.".\
                             format(module=self))
        
    def log_segment(self, segment):
        json.dump(segment, self.logfile,
                  default=self._convert_segment, 
                  separators=(',',':'),
                  sort_keys = True)
           
    def read_segment(self):
        segment = json.load(self.logfile, object_hook=self._reconvert_segment)
        return segment
        
    def print_segment(self, segment):
        print ( "segment:", segment.getNameString() )
        print ( "probability =", segment.getProbability() )
        print ( "parent =", segment.getParentNameString() )

    def _convert_bin(self, _bin):
        """
        Convert bin to built-in types
        @return list of dictionaries with bindata
        """
        self.bin_content = []
        self.dictionary = {\
            'i':_bin.getIterationId(),
            'b':_bin.getId(),
            'r':_bin.getReferenceSegmentName(),
            'n':_bin.getTargetNumberOfSegments()}
        self.bin_content.append(self.dictionary)
        for segment in _bin.segments:
            self.bin_content.append(self._convert_segment(segment))
        return self.bin_content
        
    def _reconvert_bin(self, read_dictionary):
        """
            @return returns a bin from saved data
        """
        # read bin data
        if set(('i', 'b', 'r', 'n')) <= set(read_dictionary):
            self.reference_name           = read_dictionary.get('r').split("-")
            self.newbin = Bin(\
                iteration_id              = int(read_dictionary.get('i')), 
                bin_id                    = int(read_dictionary.get('b')), 
                reference_iteration_id    = int(self.reference_name[0]), 
                reference_bin_id          = int(self.reference_name[1]),
                reference_segment_id      = int(self.reference_name[2]),
                target_number_of_segments = int(read_dictionary.get('n')) )
        
        # read segment data        
        elif set(('p', 'i', 'b', 's')) <= set(read_dictionary):
            self.newbin.generateSegment(\
                probability = read_dictionary.get('p'),
                parent_bin_id = read_dictionary.get('b'),
                parent_segment_id = read_dictionary.get('s'))            
        else:
            raise Exception("Not a bin or segment given for reconversion.")

        return self.newbin
        
    def log_bin(self, _bin):
        json.dump(_bin, self.logfile, default=self._convert_bin,
                          sort_keys=True, separators=(',',':'))
        self.logfile.write("\n")
        self.logfile.flush()
    
    def read_bin(self, line):
        _bin = json.loads(line, object_hook=self._reconvert_bin)
        #~ print("json returns:\n", _bin)
        #~ print ("newly read bin: ")
        #~ self.print_bin(_bin)
        #~ print ()
        return self.newbin
    
    def read_bins(self):
        bins = []
        for line in self.readfile:
            print (line)
            _bin = json.loads(line, object_hook=self._reconvert_bin)
            bins.append(_bin)
        return bins
        
            
    def print_bin(self, _bin):
        print ("____Bin:____")
        print ("Bin {i:05d}-{b:05d}".format(i=_bin.getIterationId(), b=_bin.getId()))
        print ("Reference: {ref}".format(ref=_bin.getReferenceSegmentName()))
        print ("Target number of segments: {t}".format(t=_bin.getTargetNumberOfSegments()))
        print ("____Segments____:")
        for segment in _bin.segments:
            self.print_segment(segment)
            
    def _convert_iteration(self):
        pass
        
    def _reconvert_iteration(self):
        pass
        
    def log_iterations(self, iterations):
        for iteration in iterations:
            for _bin in iteration.bins:
                self.log_bin(_bin)
    
    def read_iterations(self):
        self.iterations = []
        with open(self.read_filename, "r") as readfile:
            self.firstline = readfile.readline()
            self.firstbin = json.loads(self.firstline, object_hook=self._reconvert_bin)
            self.iteration_id = self.newbin.getIterationId()
            self.iteration = Iteration(self.iteration_id)
            readfile.seek(0)
            for self.line in readfile:
                """
                 for now json basically just starts self._reconvert_bin where a
                 self.newbin is created
                """
                json.loads(self.line, object_hook=self._reconvert_bin)
                # create new iteration if necessary
                if self.newbin.getIterationId() != self.iteration_id:
                    self.iterations.append(self.iteration)
                    self.iteration_id = self.newbin.getIterationId()
                    self.iteration = Iteration(self.iteration_id)
                    
                self.iteration.bins.append(self.newbin)
            self.iterations.append(self.iteration)
        return self.iterations

    def print_iteration(self, iteration):
        for _bin in iteration.bins:
            self.print_bin(_bin)
        
    def close(self):
        self.logfile.close()
        #~ if os.path.isfile(self.read_filename):
            #~ self.readfile.close()
