"""
MD module provides an interface to the langevin propagator.
"""
from __future__ import print_function
import os
import ConfigParser
import numpy
import sys
from datetime import datetime
import pickle
import config_parser
import numpy as np

# We need the hdWE program path in the PYTHONPATH when
# langevin_module is started by mpirun
if __name__ == "__main__":
    filepath  = os.path.realpath(__file__)
    dirname   = os.path.dirname(filepath)
    parentdir = os.path.dirname(dirname)
    sys.path.append(parentdir)

class MD_module():
    # WORKDIR                     = ''
    # debug                       = False
    # configuration_file          = ''   path to the amer configuration file that contains the 
    #                                    following entries:
    # langevin_topology_file      = ''   the langevin topology file
    # langevin_infile             = ''   the langevin in file
    # langevin_binary             = ''   sander, pmemd, pmemd.cuda
    # parallelization_mode        = ''   serial
    
    def __init__(self, CONFIGFILE, debug):
        """Initializes the MD module, reads the langevin configuration file
           and sets up the Langevin integrator.
        """
        self.debug                 = debug
        try:
            self.loadConfigFile(CONFIGFILE)
        except:
            os.system('sync; sleep 2')
            self.loadConfigFile(CONFIGFILE)
        
        # fire up Langevin integrator
        sys.path.append(os.path.dirname(self.langevin_class_file))
        try:
            from bbk import Langevin
            self.langevin = Langevin(self.langevin_infile)
        except:
            sys.stderr.write('Could not load Langevin integrator.\n')
            sys.exit(-1)
    
    def loadConfigFile(self, CONFIGFILE):
        #read in configuration file
        self.configfile            = CONFIGFILE
        config                     = ConfigParser.ConfigParser()
        config.read(CONFIGFILE)
        self.workdir               = config.get('hdWE','WORKDIR')
        self.jobname               = config.get('hdWE','JOBNAME')
        self.langevin_infile       = config.get('langevin','infile-path')
        self.langevin_class_file   = config_parser.parseLangevinClassFile(config)
        self.parallelization_mode  = config.get('langevin','parallelization-mode')
        self.tmp_dir               = config_parser.parseLangevinTmpDirPath(config) 
        # Only search mpirun binary in config file if MPI switched on
        if self.parallelization_mode == 'mpi':
            self.mpirun                = config.get('langevin', 'mpirun')
        
        # local link to iteration
        self.iteration              = None
        self.ITERATION_DUMP_FNAME   = '{jn}-run/iteration.dump'.format(jn=self.jobname)
        self.COORDINATES_DUMP_FNAME = '{jn}-run/coordinates.dump'.format(jn=self.jobname)
        
        self.keep_trajectory_files = bool(config.get('hdWE', 'keep-trajectory-files').lower() == "true")
        
        # check topology and infile:
        if not os.path.isfile(self.langevin_infile):
            raise Exception("No infile found at given path:\n{}".format(self.langevin_infile))
        
        langevin_config = ConfigParser.ConfigParser()
        langevin_config.read(self.langevin_infile)
        
        # Timestep for velocity calculations
        self.dt = langevin_config.getfloat('langevin', 'dt')
        self.N_DIMENSIONS         = 1
                        
    def setIteration(self, iteration):
        self.iteration = iteration
    
    def runSegmentMD(self, segment):
        """
        Function that runs one single segment MD.
        """
        langevin_start_coords_path  = "{jn}-run/{seg}.rst7".format(jn=self.jobname, seg=segment.getParentNameString())
        langevin_end_coords_path   = "{jn}-run/{seg}.rst7".format(jn=self.jobname, seg=segment.getNameString())
         
        x_coords = []
        dists    = []
        start_coords = np.loadtxt(langevin_start_coords_path)
        x_prev = start_coords[0]
        x_curr = start_coords[1]
        x_next = self.langevin.next(x_prev, x_curr)
                                
#         command_line = self.langevin_binary + \
#                                   ' -i '   + self.langevin_infile + \
#                                   ' -c '   + langevin_start_coords_path + \
#                                   ' -r '   + langevin_end_coords_path
# 
#         # Command line for debugging
#         if self.debug:
#             os.system('echo {line} >> {jn}-debug/langevin_command_lines.log'.format(line=command_line, jn=self.jobname))
                      
        # Log and Run MD
        if self.debug:
            logfile = open("{jn}-log/{it}.MD_log".format(jn=self.jobname, it=self.iteration.getNameString()),'a')
            logfile.write(self.MdLogString(segment, status = 0))
        
        # Run MD
        for i in range(self.langevin.nstlim):
            # propagate by one step
            x_prev = x_curr
            x_curr = x_next
            x_next = self.langevin.next(x_prev, x_curr)
            if i % self.langevin.nstxout == 0:
                x_coords.append(x_curr)
            if i % self.langevin.nstdout == 0:
                dists.append(self.langevin.getNorm(x_curr))
        
        # save final coordinates
        np.savetxt(langevin_end_coords_path, self.langevin.getState(x_prev, x_curr))
                
        if self.keep_trajectory_files:
            keep_trajectory_path   = "{jn}-run/{seg}.nc".format(jn=self.jobname, seg=segment.getNameString())
            np.savetxt(keep_trajectory_path, x_coords)                
                
        if self.debug:
            logfile.write(self.MdLogString(segment, status = 1 ))
            logfile.close()
        
        return 0    
    
    def MdLogString(self, segment, status):
        """Returns a string containing system time for MD run logging."""
        if status==0:
            string = 'MD run ' + segment.getNameString() + ' start: ' + str(datetime.now()) + '\n'    
        elif status==1:
            string = 'MD run ' + segment.getNameString() + ' end:   ' + str(datetime.now()) + '\n'     
        else:
            string = 'MD run ' + segment.getNameString() + ' skipped: ' + str(datetime.now()) + '\n'    
        return string 
        
    def writeMdStatus(self, segment, MD_run_count):
        """Writes the actual WE run status in a string."""
        
        active_segments = self.iteration.getNumberOfActiveSegments()
        active_bins     = self.iteration.getNumberOfActiveBins()
        bins            = self.iteration.getNumberOfBins()
        empty_bins      = self.iteration.getNumberOfEmptyBins()
        inactive_bins   = bins - active_bins
        string = '\r     Number of active (inactive/empty) bins {active_bins} ({inactive_bins}/{empty_bins}); '\
                 'Segment {segment:05d}/{active_segments:05d}'.format(
                         active_bins     = active_bins,
                         inactive_bins   = inactive_bins,
                         empty_bins      = empty_bins,
                         segment         = MD_run_count,
                         active_segments = active_segments)
        return string
    
    def printMdStatus(self, segment, MD_run_count):
        """Writes the actual WE run status to stdout."""
        sys.stdout.write(self.writeMdStatus(segment, MD_run_count))
        sys.stdout.flush()  
    
    def runMDs(self, iteration):
        """Propagates the trajectories corresponding to an iteration using langevin."""
        self.iteration = iteration
        #Serial Run
        if self.parallelization_mode=='serial':
            MD_run_count  = 0
            for bin_loop in self.iteration:
                for segment_loop in bin_loop:
                    MD_run_count  += 1
                    self.runSegmentMD(segment_loop)
                    self.printMdStatus(segment_loop, MD_run_count)
        else:
            sys.stderr.flush()
            sys.stderr.write('Error: No parallelization implemented. Use serial mode.')
            sys.exit(-1)
    
    def distance(self, x):
        return numpy.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    
    def calcSegmentCoordinates(self, segment):
        """
        Calculates the coordinates of a segment with respect to defined dimensions
        sets segment.coordinates to calculated ones
        @return list of float coordinates 
        """
        segment_name_string = segment.getNameString()
        rst7_file = "{jn}-run/{namestring}.rst7".format(jn=self.jobname, namestring = segment_name_string)               
        coordinates = self.calcCoordinatesOfFile(rst7_file)
        segment.setCoordinates(coordinates)
        return coordinates
    
    def calcCoordinatesOfFile(self, filename):
        """
        Takes a file with optional relative path (like the run directory).
        Calculates the coordinates of a file with respect to defined dimensions.
        @return list of float coordinates 
        """
        x = numpy.loadtxt(filename)
        coordinates = [self.distance(x[1])]
        return coordinates
        
    
    def calcSegmentVelocities(self, segment):
        
        segment_name_string = segment.getNameString()
        
        rst7_file = "{jn}-run/{namestring}.rst7".format(jn=self.jobname, namestring = segment_name_string)               
        x = numpy.loadtxt(rst7_file)
        
        dx = x[1] - x[0]
        v  = dx/self.dt
        
        segment.setVelocities(v)
        return v
   
    def loadIterationFromDumpFile(self, retry=False):
        """
        Load the local copy of iteration from the 
        default dump file name 
        """
        try:
            iteration_dump_file = open(self.ITERATION_DUMP_FNAME, "r")
            self.setIteration(pickle.load(iteration_dump_file))
            iteration_dump_file.close()
        except:
            # Sync, wait, retry
            if retry == False:
                os.system('sync; sleep 2;')
                self.loadIterationFromDumpFile(retry=True)
            else:
                raise
    
    def dumpIterationToFile(self):
        """
        Dump the local iteration object to file
        """
        iteration_dump_file = open(self.ITERATION_DUMP_FNAME, 'w')
        pickle.dump(self.iteration, iteration_dump_file, protocol=2)
        iteration_dump_file.close()
    
    def removeIterationDumpFile(self):
        """
        Delete the default iteration dump file
        """
        os.remove(self.ITERATION_DUMP_FNAME)
    
    def dumpCoordinatesToFile(self, coordinate_ids):
        """
        Store the coordinate matrix to dump file
        """        
        coordinate_ids_dump_file = open(self.COORDINATES_DUMP_FNAME , 'w')
        pickle.dump(coordinate_ids, coordinate_ids_dump_file)
        coordinate_ids_dump_file.close()

    def loadCoordinatesFromDumpFile(self):
        """
        Load the coordinate matrix from dump file
        """
        coordinates_dump_file = open(self.COORDINATES_DUMP_FNAME , 'r')
        coordinates = pickle.load(coordinates_dump_file)
        coordinates_dump_file.close()
        return coordinates

    def removeCoordinatesDumpFile(self):
        """ Remove coordinate matrix dump file """
        os.remove(self.COORDINATES_DUMP_FNAME)
    

    def calcCoordinates(self, iteration):
        """
        Calculates the coordinates with respect 
        to binned dimensions for all segments in iteration
        sets these coordinates in each Segment
        @return matrix of coordinates [N_segments x N_dimensions]
        """
        self.setIteration(iteration)
        if self.parallelization_mode == "serial":
            # calculate and set coordinates
            coordinates = numpy.zeros([iteration.getNumberOfSegments(), self.N_DIMENSIONS])
            i = 0
            for this_bin in iteration:
                for this_segment in this_bin:
                    self.calcSegmentCoordinates(this_segment)
                    for j in range(self.N_DIMENSIONS):
                        coordinates[i,j] = this_segment.getCoordinates()[j]
                    i += 1
        else:
            sys.stderr.flush()
            sys.stderr.write('Error: No parallel version available. Use serial mode.')
            sys.exit(-1)
        
        return coordinates

    def calcVelocity(self, iteration):
        """
        Calculate the relative velocities (at the moment only 
        the 2 ions special case is supported)
        """
        self.setIteration(iteration)
        if self.parallelization_mode == "serial":
            # calculate and set coordinates
            for this_bin in iteration:
                for this_segment in this_bin:
                    self.calcSegmentVelocities(this_segment)
        
        else:
            sys.stderr.flush()
            sys.stderr.write("ERROR: Velocity calculation is only available in serial mode")
            sys.exit(-1)
            
    def removeCoordinateFiles(self, iteration, compress, keep=0):
        """
        @param keep: do not delete the first #keep segment coordinate files per bin
        """
        for this_bin in iteration:
            for this_segment in this_bin.segments[keep:]:
                    try:
                        os.remove("{jn}-run/{seg}.rst7".format(jn=self.jobname, seg=this_segment.getNameString()))
                    except OSError:
                        continue
                