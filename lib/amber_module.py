"""
MD module provides an interface to the AMBER molecular dynamics
toolchain.
"""
from __future__ import print_function
import os
import ConfigParser
import numpy
import sys
from datetime import datetime
import pickle
import uuid
import config_parser

# We need the hdWE program path in the PYTHONPATH when
# amber_module is started by mpirun
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
    # amber_topology_file         = ''   the amber topology file
    # amber_infile                = ''   the amber in file
    # amber_coordinate_mask       = ''   for example for RMSD calculation. Example ':1-3@CA'
    # amber_binary                = ''   sander, pmemd, pmemd.cuda
    # parallelization_mode        = ''   serial, mpi
    
    def __init__(self, CONFIGFILE, debug):
        """Initializes the MD module and Reads the amber configuration file.
        """
        self.debug                 = debug
        try:
            self.loadConfigFile(CONFIGFILE)
        except:
            os.system('sync; sleep 2')
            self.loadConfigFile(CONFIGFILE)
    
    def loadConfigFile(self, CONFIGFILE):
        #read in configuration file
        self.configfile            = CONFIGFILE
        config                     = ConfigParser.ConfigParser()
        config.read(CONFIGFILE)
        self.workdir               = config.get('hdWE','WORKDIR')
        self.jobname               = config.get('hdWE','JOBNAME')
        self.amber_topology_file   = config.get('amber','topology-path')
        self.amber_infile          = config.get('amber','infile-path')
        self.amber_binary          = config_parser.parseAmberBinary(config)
        self.cpptraj_binary        = config_parser.parseCpptrajBinary(config)
        self.coordinate_masks_file = config.get('amber', 'coordinate-masks')
        self.rmsd_mask             = config_parser.parseAmberRmsdMask(config)
        self.rmsd_fit_mask         = config_parser.parseAmberRmsdFitMask(config)
        self.parallelization_mode  = config.get('amber','parallelization-mode')
        # Only search mpirun binary in config file if MPI switched on
        if self.parallelization_mode == 'mpi':
            self.mpirun                = config.get('amber', 'mpirun')
        
        # Check for GPU support
        self.has_cuda = (os.path.basename(self.amber_binary) == "pmemd.cuda")
        if self.has_cuda: 
            # A list of cuda visible devices indices [0,1,...]
            self.cuda_visible_devices = map(int, config.get('amber', 'cuda_visible_devices').strip().split(','))
        # local link to iteration
        self.iteration              = None
        self.ITERATION_DUMP_FNAME   = '{jn}-run/iteration.dump'.format(jn=self.jobname)
        self.COORDINATES_DUMP_FNAME = '{jn}-run/coordinates.dump'.format(jn=self.jobname)
        self.BINS_RMSDS_DUMP_FNAME  = '{jn}-run/binsrmsds.dump'.format(jn=self.jobname)
        
        self.keep_trajectory_files = bool(config.get('hdWE', 'keep-trajectory-files').lower() == "true")
        
        # check topology and infile:
        if not os.path.isfile(self.amber_topology_file):
            raise Exception("No topology found at given path:\n{}".format(self.amber_topology_file))
        if not os.path.isfile(self.amber_infile):
            raise Exception("No infile found at given path:\n{}".format(self.amber_infile))
        
        # get coordinate masks
        self.COORDINATE_MASKS = []
        c_mask_file = open(self.coordinate_masks_file,'r')
        for line in c_mask_file.readlines():
            self.COORDINATE_MASKS.append(line.strip())
        self.N_DIMENSIONS = len(self.COORDINATE_MASKS)
        
    def setIteration(self, iteration):
        self.iteration = iteration
    
    def runSegmentMD(self, segment):
        """
        Function that runs one single segment MD.
        """

        amber_start_coords_path = "{jn}-run/{seg}.rst7".format(jn=self.jobname, seg=segment.getParentNameString()) 
        amber_end_coords_path   = "{jn}-run/{seg}.rst7".format(jn=self.jobname, seg=segment.getNameString()) 
        
        if self.debug:
            amber_info_path         = "{jn}-run/{seg}.inf".format(jn=self.jobname, seg=segment.getNameString())
            amber_outfile_path      = "{jn}-run/{seg}.out".format(jn=self.jobname, seg=segment.getNameString())
            amber_trajectory_path   = "{jn}-run/{seg}.nc".format( jn=self.jobname, seg=segment.getNameString())
        else:
            #amber_info_path         = '/dev/null'
            #amber_trajectory_path   = '/dev/null'
            #amber_outfile_path      = '/dev/null'
            #TODO: Strange bug on REX allows only info and traj to be /dev/null
            amber_info_path         = '/tmp/{seg}_{id}.inf'.format(seg=segment.getNameString(), id=uuid.uuid1())
            amber_outfile_path      = '/tmp/{seg}_{id}.out'.format(seg=segment.getNameString(), id=uuid.uuid1()) 
            amber_trajectory_path   = '/tmp/{seg}_{id}.nc'.format( seg=segment.getNameString(), id=uuid.uuid1()) 

        #TODO dou ble gpu usage on one node, can maybe by done more elegantly?

        if self.keep_trajectory_files:
            amber_trajectory_path   = "{jn}-run/{seg}.nc".format(jn=self.jobname, seg=segment.getNameString())
        
        if self.has_cuda:
            # Assign one cuda visible device per job
            os.environ['CUDA_VISIBLE_DEVICES'] = str(self.cuda_visible_devices[segment.getId() % len(self.cuda_visible_devices)])
 
        # Overwrite the previous settings                       
        command_line = self.amber_binary + ' -O' + \
                                  ' -p '   + self.amber_topology_file + \
                                  ' -i '   + self.amber_infile + \
                                  ' -c '   + amber_start_coords_path + \
                                  ' -o '   + amber_outfile_path + \
                                  ' -x '   + amber_trajectory_path + \
                                  ' -inf ' + amber_info_path + \
                                  ' -r '   + amber_end_coords_path
                                    
        

        #Command line for debugging
        if self.debug:
            os.system('echo {line} >> {jn}-debug/amber_command_lines.log'.format(line=command_line, jn=self.jobname))
                      
        #Log and Run MD
        if self.debug:
            logfile = open("{jn}-log/{it}.MD_log".format(jn=self.jobname, it=self.iteration.getNameString()),'a')
            logfile.write(self.MdLogString(segment, status = 0))
        
        # Run MD
        ret_code = os.system(command_line)
        
        # If anything went wrong, sync, wait a second, then retry
        if ret_code != 0:
            os.system('sync; sleep 1')
            ret_code = os.system(command_line)
        
        if self.debug:
            logfile.write(self.MdLogString(segment, status = 1 ))
            logfile.close()
        
        # Remove out files
        if not self.debug:
            try: os.remove(amber_outfile_path)
            except OSError: pass
            try: os.remove(amber_info_path)
            except OSError: pass
            if not self.keep_trajectory_files:
                try: os.remove(amber_trajectory_path)
                except OSError: pass
    
    def SkipSegmentMD(self, segment, MD_run_count, MD_skip_count):
        """Function that runs one single segment MD."""
        command_line = self.skipCommandLineString(segment)
        #Command line for debugging
        if self.debug:
                os.system('echo ' + command_line + \
                          ' >> {jn}-debug/amber_command_lines.log'.format(jn=self.jobname))
        #Log and Run MD
        logfile = open("{jn}-log/{it}.MD_log".format(jn = self.jobname,
                                                        it = self.iteration.getNameString()),
                       'a')
        logfile.write(self.MdLogString(segment, status = 2 ))
        #sys.stdout.write(self.writeMdStatus(segment, MD_run_count, MD_skip_count))
        #sys.stdout.flush()
        os.system(command_line)
        logfile.close()
    
    def skipCommandLineString(self, segment):
        """
        Returns the command line for linking segment restart files of skipped bins to next iteration.
        """
        amber_start_coords_path = "{jn}-run/{pseg}.rst7".format(jn=self.jobname, pseg=segment.getParentNameString())
        amber_end_coords_path = "{jn}-run/{seg}.rst7".format(jn=self.jobname, seg=segment.getNameString())        
        
        skip_command_line       = 'ln {start} {end}'.format(
                                                  start=amber_start_coords_path,
                                                  end=amber_end_coords_path)
        return skip_command_line
      
    def MdLogString(self, segment,status):
        """Returns a string containing system time for MD run logging."""
        if status==0:
            string = 'MD run ' + segment.getNameString() + ' start: ' + str(datetime.now()) + '\n'    
        elif status==1:
            string = 'MD run ' + segment.getNameString() + ' end:   ' + str(datetime.now()) + '\n'     
        else:
            string = 'MD run ' + segment.getNameString() + ' skipped: ' + str(datetime.now()) + '\n'    
        return string 
        
    def writeMdStatus(self, segment, MD_run_count, MD_skip_count):
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
    
    def printMdStatus(self, segment, MD_run_count, MD_skip_count):
        """Writes the actual WE run status to stdout."""
        sys.stdout.write(self.writeMdStatus(segment, MD_run_count, MD_skip_count))
        sys.stdout.flush()  
    
    def runMDs(self, iteration):
        """Propagates the trajectories corresponding to an iteration using amber."""
        self.iteration = iteration
        #Serial Run
        if self.parallelization_mode=='serial':
            MD_run_count  = 0
            MD_skip_count = 0
            for bin_loop in self.iteration:
                for segment_loop in bin_loop:
                    MD_run_count  += 1
                    self.runSegmentMD(segment_loop)
                    self.printMdStatus(segment_loop, MD_run_count, MD_skip_count)
                                      
        #MPI parallel run
        if self.parallelization_mode == 'mpi':
            # dump the iteration object to a file
            self.dumpIterationToFile()
            # Call myself with MPI...
            if self.debug:
                print("Switching to MPI for MD now...")
            error_code = os.system("{0} {1} {2} {3} {4} {5} ".format(self.mpirun, 
                                  sys.executable,
                                  os.path.realpath(__file__),
                                  "MPIMD",
                                  self.configfile,
                                  str(self.debug)))
            if error_code != 0:
                sys.stderr.flush()
                sys.stderr.write("Received error code during runMDs: {0}\n".format(error_code))
                sys.exit()

    def calcSegmentCoordinates(self, segment):
        """
        Calculates the coordinates of a segment with respect to defined dimensions
        sets segment.coordinates to calculated ones
        @return list of float coordinates 
        """
        # Write the cpptraj infile
        segment_name_string = segment.getNameString()
       
        coordinates = self.calcCoordinatesOfFile("{jn}-run/{namestring}.rst7".format(jn=self.jobname, 
                                                                                     namestring = segment_name_string))
        # set coordinates in segment
        segment.setCoordinates(coordinates)
        return coordinates

    def calcCoordinatesOfFile(self, filename):
        """
        Takes a file with optional relative path (like the run directory).
        Calculates the coordinates of a file with respect to defined dimensions.
        Sets segment.coordinates to calculated ones.
        @return list of float coordinates 
        """
        # Write the cpptraj infile
        UUID = uuid.uuid1()
        if self.debug:
            cpptraj_infile_path  = "{filename}.cpptraj_in".format(jn=self.jobname, filename = filename)
            cpptraj_outfile_path = "{filename}.cpptraj_out".format(jn=self.jobname, filename = filename)
            cpptraj_logfile_path = "{jn}-log/cpptraj.log".format(jn=self.jobname)

        else:
            basename = filename.split('/')[-1]
            cpptraj_infile_path = "/tmp/{0}_{1}.cpptraj_in".format(basename, UUID)
            cpptraj_outfile_path = "/tmp/{0}_{1}.cpptraj_out".format(basename, UUID)
            cpptraj_logfile_path = "/dev/null"

        
        cpptraj_infile      = open(cpptraj_infile_path, 'w')
        cpptraj_infile.write('trajin {filename}\n'.format(filename = filename))       
        for coordinate_mask in self.COORDINATE_MASKS:
                cpptraj_infile.write('{mask} out {out}\n'.format(mask = coordinate_mask,
                                                                 out  = cpptraj_outfile_path))
        
        # Write changes to file
        cpptraj_infile.close()
        
        # Run cpptraj
        cpptraj_execute_string = ' -p {top} -i {inpath} > {log}'.format(
                                                            top    = self.amber_topology_file, 
                                                            inpath = cpptraj_infile_path,
                                                            log     = cpptraj_logfile_path)
        os.system('{cpptraj} {execute}'.format(cpptraj=self.cpptraj_binary, execute=cpptraj_execute_string))
                
        # Load cpptraj output as numpy array
        try:
            coordinates = numpy.loadtxt(cpptraj_outfile_path) 
            # Delete the first entry which refers to the frame index
            coordinates = numpy.delete(coordinates, 0)
        except:
            #TODO What should happen then?
            sys.stderr.write('amber_module error: cpptraj output {0} can not '\
                  'be found or loaded.\n'.format(cpptraj_outfile_path))
        
        if not self.debug:
            os.remove(cpptraj_outfile_path)
            os.remove(cpptraj_infile_path)
        
        return coordinates
    
    def getRmsdsToSegments(self, segments, target_number_of_segments=0):
        """
        Calculates the rmsd of all segments with respect to each other.
        @return matrix of rmsds
        """
        # Only work if merge is required
        if len(segments) <= target_number_of_segments:
            return 0.0
        # Use ParentNameString because this routine runs before the MD
        
        # Write the cpptraj infile
        segment_name_string = segments[0].getParentNameString() 
        UUID = uuid.uuid1()
        
        if self.debug:
            cpptraj_infile_path  = "{jn}-run/{segment}.cpptraj_rmsds_in".format(jn=self.jobname, segment=segment_name_string)
            cpptraj_outfile_path = "{jn}-run/{segment}.cpptraj_rmsds_out".format(jn=self.jobname,segment=segment_name_string)
            cpptraj_logfile_path = "{jn}-log/cpptraj_rmsds.log".format(jn=self.jobname)

        else:
            cpptraj_infile_path  = "/tmp/{0}_{1}.cpptraj_rmds_in".format(segment_name_string, UUID)
            cpptraj_outfile_path = "/tmp/{0}_{1}.cpptraj_rmds_out".format(segment_name_string, UUID)
            cpptraj_logfile_path = "/dev/null"

        
        # cpptraj input
        cpptraj_infile = open(cpptraj_infile_path, 'w')
        for this_segment in segments:
            cpptraj_infile.write('trajin {jn}-run/{segment}.rst7\n'.format(jn      = self.jobname, 
                                                                           segment = this_segment.getParentNameString()))                   
            cpptraj_infile.write('reference {jn}-run/{ref}.rst7 [{refname}]\n'.format(jn=self.jobname,
                                                                    ref     = this_segment.getParentNameString(),
                                                                    refname = this_segment.getParentNameString() ))
            # Aligning
            cpptraj_infile.write('rmsd {refname}.fit {fit_mask} ref [{refname}]\n'.format(fit_mask = self.rmsd_fit_mask,
                                                                                 out      = cpptraj_outfile_path,
                                                                                 refname  = this_segment.getParentNameString()))
            # Actual RMSD calculation
            cpptraj_infile.write('rmsd {refname} {mask} nofit out {out} ref [{refname}]\n'.format(mask = self.rmsd_mask,
                                                                                 out      = cpptraj_outfile_path,
                                                                                 refname  = this_segment.getParentNameString()))
        
        cpptraj_infile.close()
        
        # Run cpptraj
        cpptraj_execute_string = ' -p {top} -i {inpath} > {log}'.format(
                                                            top    = self.amber_topology_file, 
                                                            inpath = cpptraj_infile_path,
                                                            log     = cpptraj_logfile_path)
        os.system('{cpptraj} {execute}'.format(cpptraj=self.cpptraj_binary, execute=cpptraj_execute_string))
                
        # Load cpptraj output as numpy array
        try:
            rmsds = numpy.loadtxt(cpptraj_outfile_path) 
            # Delete the first entry which refers to the frame index
            rmsds = rmsds[:,1:]
        except:
            #TODO What should happen then?
            sys.stderr.write('amber_module error: cpptraj output {0} can not '\
                  'be found or loaded.\n'.format(cpptraj_outfile_path))
        
        if not self.debug:
            os.remove(cpptraj_outfile_path)
            os.remove(cpptraj_infile_path)

        return rmsds

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


    def dumpBinsRmsdsToFile(self, bins_rmsds):
        """
        Store the bins segment vs. segments rmsds array from dump file
        """
        bins_rmsds_dump_file = open(self.BINS_RMSDS_DUMP_FNAME, 'w')
        pickle.dump(bins_rmsds, bins_rmsds_dump_file)
        bins_rmsds_dump_file.close()
        
    def loadBinsRmsdsFromDumpFile(self):
        """
        Load the bins segment vs. segments rmsds array from dump file
        """
        bins_rmsds_dump_file = open(self.BINS_RMSDS_DUMP_FNAME , 'r')
        bins_rmsds = pickle.load(bins_rmsds_dump_file)
        bins_rmsds_dump_file.close()
        return bins_rmsds
    
    def removeBinsRmsdsToFile(self):
        """ Remove the bins segment vs. segments rmsds array from dump file """
        os.remove(self.BINS_RMSDS_DUMP_FNAME)

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
                    # write to coordinates matrix 
                    # (not necessary because calcSegmentCoordinates sets them in the segment)
                    for j in range(self.N_DIMENSIONS):
                        coordinates[i,j] = this_segment.getCoordinates()[j]
                    i += 1
        
        if self.parallelization_mode == "mpi":
            # dump the iteration object to a file
            self.dumpIterationToFile()
            # Call myself with MPI...
            if self.debug:
                print("Switching to MPI for coordinate calculation.")
            error_code = os.system("{0} {1} {2} {3} {4} {5} ".format(self.mpirun, 
                                                        sys.executable,
                                                        os.path.realpath(__file__),
                                                        "MPICOORD",
                                                        self.configfile, 
                                                        str(self.debug)))
            if error_code != 0:
                sys.stderr.flush()
                sys.stderr.write("Received error code during calcCoordinates: {0}\n".format(error_code))
                sys.exit()
            coordinates = self.loadCoordinatesFromDumpFile()
            if not self.debug:
                self.removeCoordinatesDumpFile()
            # assign coordinates to segments
            i=0
            for this_bin in iteration:
                for this_segment in this_bin:
                    this_segment.setCoordinates(coordinates[i])
                    i += 1
        return coordinates

    def calcBinsRmsds(self, iteration):
        """
        Calculates the segments vs segments rmsds for each  
        bin.
        @return array of matrices with segment vs. segment rmsds [N_bins x [N_segments x N_segments]]
        """
        self.setIteration(iteration)
        if self.parallelization_mode == "serial":
            # Setup the list of matrices
            rmsds = numpy.array([None] * iteration.getNumberOfBins(), dtype=object)
            i = 0
            for this_bin in iteration:
                rmsds[i] = self.getRmsdsToSegments(this_bin.segments, this_bin.getTargetNumberOfSegments())
                i += 1
        
        if self.parallelization_mode == "mpi":
            # dump the iteration object to a file
            self.dumpIterationToFile()
            # Call myself with MPI...
            if self.debug:
                print("Switching to MPI for bins rmsd matrix calculation.")
            error_code = os.system("{0} {1} {2} {3} {4} {5} ".format(self.mpirun, 
                                                        sys.executable,
                                                        os.path.realpath(__file__),
                                                        "MPIBINRMSDS",
                                                        self.configfile, 
                                                        str(self.debug)))
            if error_code != 0:
                sys.stderr.flush()
                sys.stderr.write("Received error code during calcBinsRmsds: {0}\n".format(error_code))
                sys.exit()
            rmsds = self.loadBinsRmsdsFromDumpFile()
            if not self.debug:
                self.removeBinsRmsdsToFile()
            # assign coordinates to segments
                
        return rmsds
        
    def ana_calcCoordinateOfSegment(self, segment_name_string, cpptraj_lines, use_trajectory):
        """
        Calculates the value of a coordinate corresponding to a segment and defined
        in cpptraj_line via cpptraj.
        """
        
        #Write the cpptraj infile
        cpptraj_infile_path = "{seg}.ana_calculatePMF_cpptraj_in".format(seg=segment_name_string)
        cpptraj_output_path = "{seg}.ana_calculatePMF_cpptraj_output".format(seg=segment_name_string)
        cpptraj_infile      = open(cpptraj_infile_path,'w')
        if use_trajectory == False:
            cpptraj_infile.write('trajin {jn}-run/{segment}.rst7\n'.format(jn=self.jobname,
                                                                           segment=segment_name_string) )
        else:
            cpptraj_infile.write('trajin {jn}-run/{segment}.nc\n'.format(jn=self.jobname,
                                                                         segment=segment_name_string) )            
        cpptraj_infile.writelines(cpptraj_lines + ' out ' + cpptraj_output_path )
        cpptraj_infile.close()
        
        #Execute cpptraj
        cpptraj_execute_string =' -p ' + self.amber_topology_file + \
                                ' -i ' + cpptraj_infile_path + \
                                ' >> {jn}-log/ana_calculatePMF_cpptraj.log'.format(jn=self.jobname)
              
        os.system('{cpptraj} {execute}'.format(cpptraj=self.cpptraj_binary, execute=cpptraj_execute_string))
        #Load cpptraj output as numpy array
        try:
            coordinates = numpy.loadtxt(cpptraj_output_path) 
        except:
            #TODO What should happen then?
            print('amber_module error: cpptraj output ' + cpptraj_output_path + ' can not be found or loaded.')
            
        #Remove temporary files
        if not self.debug:
            os.remove(cpptraj_infile_path)
            os.remove(cpptraj_output_path)
         
        if use_trajectory == False:   
            return  [coordinates[1]]
        else:
            return coordinates[:,1]
    
    def removeCoordinateFiles(self, iteration, compress, keep=0):
        """
        @param keep: do not delete the first #keep segment coordinate files per bin
        """
        # to avoid data loss, check if the compressed .nc file has been created.
        # do not remove the rst7 files in case an error occured during compression  
        if compress == True:
            nc_file_path = '{jn}-run/{iterationId}.nc\n'.format(jn = self.jobname, iterationId=iteration.getNameString())
            if os.path.isfile(nc_file_path) == False:
                    print('No compressed trajectory file (.nc) of iteration {:05d} found although compression was requested.'.format(iteration.getNameString())) 
                    print('Not deleting .rst7 files.')
                    return

        for this_bin in iteration:
            for this_segment in this_bin.segments[keep:]:
                    try:
                        os.remove("{jn}-run/{seg}.rst7".format(jn=self.jobname, seg=this_segment.getNameString()))
                    except OSError:
                        continue
                

    def compressIteration(self, iteration, cpptraj_closest_mask=""):
        """
        puts all .rst7 files from an iteration together into a .nc binary file, omitting velocities.
        a segment_list is created from which the original rst7 files can be recreated.
        """

        # get all segment name strings in a list
        segment_name_list = []
        for this_bin in iteration:
            for this_segment in this_bin:
                segment_name_list.append(this_segment.getNameString())
        
        # write the segment names into an index file
        segment_name_list_path = "{jn}-run/{iterationId}.segment_list".format(jn=self.jobname,
                                                                              iterationId=iteration.getNameString())
        segment_name_list_file = open(segment_name_list_path, 'w')
        for segment_name in segment_name_list:
            segment_name_list_file.write(segment_name + '\n')
        segment_name_list_file.close()
        
        # write cppraj in file
        UUID                 = uuid.uuid1()
        cpptraj_file_path    = "/tmp/{iterationId}_{uuid}.compress_cpptraj_in".format(iterationId=iteration.getNameString(), uuid=UUID)
        cpptraj_logfile_path = "/dev/null"
        cpptraj_file         = open(cpptraj_file_path, 'w')    
        for segment_name in segment_name_list:
            cpptraj_file.write('trajin {jn}-run/{segment_name}.rst7\n'.format(jn = self.jobname, segment_name = segment_name))
        if cpptraj_closest_mask != "":
            cpptraj_file.write('closest {mask}\n'.format(mask=cpptraj_closest_mask))
            cpptraj_file.write('autoimage\n')
        cpptraj_file.write('trajout {jn}-run/{iterationId}.nc netcdf novelocity\n'.format(jn = self.jobname, iterationId=iteration.getNameString()))
        cpptraj_file.close()
        
        # execute cpptraj
        cpptraj_execute_string = ' -p {top} -i {inpath} > {log}'.format(
                                                            top    = self.amber_topology_file, 
                                                            inpath = cpptraj_file_path,
                                                            log     = cpptraj_logfile_path)
        os.system('{cpptraj} {execute}'.format(cpptraj=self.cpptraj_binary, execute=cpptraj_execute_string))
        
        # remove cpptraj infile
        os.remove(cpptraj_file_path)
        
    
    
###########################################
#     WHEN CALLED IN MPI CONTEXT          #
###########################################
def doMPIMD(CONFIGFILE, debug):
    """
    This function can run on multiple MPI processes in parallel
    to propagate all segments of iteration via MD 
    """    
    # Read in the call arguments
    CONFIGFILE               = CONFIGFILE
    debug                    = bool(debug == "True")
    # Initialize MD module
    md_module = MD_module(CONFIGFILE = CONFIGFILE, debug = debug)
    md_module.loadIterationFromDumpFile()
    if debug:
        if rank == 0:
            sys.stderr.write("Number of MPI processes: {0}\n".format(size))
            sys.stderr.flush()
        comm.barrier()
    # Every node works only on segments modulo their rank 
    workcount = 0
    md_skip_count = 0
    for loop_bin in md_module.iteration:
        for loop_segment in loop_bin:
            #if not loop_bin.isConverged():
            if workcount % size == rank:
                # Run MD on this node
                md_module.runSegmentMD(loop_segment)
            #else:
            #    md_skip_count += 1
            #    if workcount % size == rank:
            #        # Run MD skip
            #        md_module.SkipSegmentMD(loop_segment, workcount, md_skip_count)
            workcount += 1
            # Log if rank 0
            if rank == 0:
                md_module.printMdStatus(loop_segment, workcount, md_skip_count)
                
    # Wait for all processes to finish
    comm.barrier()
    if rank == 0:
        md_module.printMdStatus(loop_segment, workcount, md_skip_count)
        #sys.stdout.write("\n")
        if debug:
            sys.stdout.write("Finishing MPI\n")
        sys.stdout.flush()
        # Remove the iteration dump file
        if not debug:
            md_module.removeIterationDumpFile()

def doMPICalcCoordinates(CONFIGFILE, debug):
    # Read in the call arguments
    CONFIGFILE               = CONFIGFILE
    debug                    = bool(debug == "True")
    # Initialize MD module
    md_module = MD_module(CONFIGFILE = CONFIGFILE, debug = debug)
    md_module.loadIterationFromDumpFile()
    iteration = md_module.iteration
    
    if debug and rank == 0:
        sys.stderr.write("Number of MPI processes: {0}\n".format(size))
        sys.stderr.flush()

    # Setup the matrix of coordinates
    coordinates = numpy.zeros([iteration.getNumberOfSegments(), md_module.N_DIMENSIONS], dtype = float)
    # calculate entries
    i = 0
    for this_bin in iteration:
        for this_segment in this_bin:
            if i % size == rank:
                segment_coordinates = md_module.calcSegmentCoordinates(this_segment)
                # fill matrix
                for j in range(md_module.N_DIMENSIONS):
                    coordinates[i][j] = segment_coordinates[j]
            i += 1
    # Collect the matrix on the root process (rank == 0)
    coordinates = comm.reduce(coordinates, op = MPI.SUM, root=0)
    # Dump matrix to file
    if rank == 0:
        md_module.dumpCoordinatesToFile(coordinates)
        if debug:
            print("All Coordinate Ids: \n", coordinates)
            print("Finished MPI Coordinates calculation")

def doMPICalcBinsRmsds(CONFIGFILE, debug):
    # Read in the call arguments
    CONFIGFILE = CONFIGFILE
    debug      = bool(debug == "True")
    # Initialize MD module
    md_module  = MD_module(CONFIGFILE = CONFIGFILE, debug = debug)
    md_module.loadIterationFromDumpFile()
    iteration  = md_module.iteration
    
    if debug and rank == 0:
        sys.stderr.write("Number of MPI processes: {0}\n".format(size))
        sys.stderr.flush()

    # Setup the list of matrices
    bins_rmsds = numpy.array([0.0] * iteration.getNumberOfBins(), dtype=object)
    
    # Select only those bins which require merge
    merge_requiring_bins = []
    for this_bin in iteration:
        if this_bin.getNumberOfSegments() > this_bin.getTargetNumberOfSegments():
            merge_requiring_bins.append(this_bin)
    # calculate entries
    i = 0
    for this_bin in merge_requiring_bins:
        if i % size == rank:
            bins_rmsds[this_bin.getId()] = md_module.getRmsdsToSegments(this_bin.segments, this_bin.getTargetNumberOfSegments())
        i += 1
        
    # Collect the matrix on the root process (rank == 0)
    bins_rmsds = comm.reduce(bins_rmsds, op = MPI.SUM, root = 0)
    # Dump matrix to file
    if rank == 0:
        md_module.dumpBinsRmsdsToFile(bins_rmsds)
        if debug:
            print("All bins RMSD matrices: \n", bins_rmsds)
            print("Finished MPI bin RMSD matrix calculation")

class MD_analysis_module():
    """ a MD module for analysis operations
        does not need the CONFIGFILE
    """
    def __init__(self, WORKDIR, JOBNAME, debug):
        """
        Initializes the MD analysis module without a config file.
        Additional input like topology or a structure might be needed for further functions.
        """
        #read in parameters
        self.workdir               = WORKDIR
        self.jobname               = JOBNAME

    
    def isSegmentFile(self, segment):
        """
        Returns true if the segment file exists in the run folder on the filesystem
        """
        return os.path.isfile("{wd}/{jn}-run/{seg}.rst7".format( wd=self.workdir, jn=self.jobname, seg=segment.getNameString()))
        
# Call only in MPI mode
if __name__ == "__main__":
    from mpi4py import MPI    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    if sys.argv[1] == "MPIMD":
        doMPIMD(CONFIGFILE=sys.argv[2], 
                debug=sys.argv[3])
    if sys.argv[1] == "MPICOORD":
        doMPICalcCoordinates(CONFIGFILE=sys.argv[2],
                             debug=sys.argv[3])
    if sys.argv[1] == "MPIBINRMSDS":
        doMPICalcBinsRmsds(CONFIGFILE=sys.argv[2],
                           debug=sys.argv[3])
