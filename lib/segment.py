#
# This file is part of hdWE. 
# Copyright (C) 2016 Manuel Luitz <manuel.luitz@tum.de>
# Copyright (C) 2016 Rainer Bomblies <r.bomblies@tum.de>
# Copyright (C) 2016 Fabian Zeller
#
# hdWE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# hdWE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with hdWE. If not, see <http://www.gnu.org/licenses/>.
# 
class Segment(object):
    """
    Defines a trajectory element of a bin
    """
    def __init__(self, probability, parent_iteration_id, parent_bin_id, parent_segment_id,
                 iteration_id, bin_id, segment_id):
        self.probability         = probability          # float
        self.parent_iteration_id = parent_iteration_id  # int
        self.parent_bin_id       = parent_bin_id        # int
        self.parent_segment_id   = parent_segment_id    # int
        self.bin_id              = bin_id               # int
        self.segment_id          = segment_id           # int
        self.iteration_id        = iteration_id         # int
        self.coordinates         = None                 # list of floats
        self.velocities          = None                 # list of floats

    def __getNameString(self, iteration_id, bin_id, segment_id):
        """
        @return Formated segment file name string
        """
        return "{iteration:08d}_{_bin:05d}_{segment:05d}"\
                .format(iteration=iteration_id,
                       _bin=bin_id,
                       segment=segment_id)
    
    def setSegmentId(self, segment_id):
        """
        required for resampling function to reorder the segment ids. 
        Should not be used from outside Bin class
        """
        self.segment_id = segment_id
        
    def getNameString(self):
        """
        @return the indices in a string following the scheme iteration_bin_segment
        """
        return self.__getNameString(iteration_id=self.iteration_id,
                                    bin_id=self.bin_id,
                                    segment_id=self.segment_id)
        
    def getParentNameString(self):
        """
        @return the indices in a string following the scheme iteration_bin_segment
        """
        return self.__getNameString(iteration_id=self.parent_iteration_id,
                                    bin_id=self.parent_bin_id,
                                    segment_id=self.parent_segment_id)
        
    def getProbability(self):
        return self.probability
        
    def setProbability(self, probability):
        self.probability = probability

    def addProbability(self, probability):
        self.probability += probability

    def subProbability(self, probability):
        self.probability -= probability

    def multProbability(self, factor):
        self.probability *= factor
        
    def getId(self):
        return self.segment_id
    
    def getBinId(self):
        return self.bin_id

    def getIterationId(self):
        return self.iteration_id

    def getParentIterationId(self):
        return self.parent_iteration_id   
                
    def getParentBinId(self):
        return self.parent_bin_id
        
    def getParentSegmentId(self):
        return self.parent_segment_id
    
    def getCoordinates(self):
        return self.coordinates
    
    def setCoordinates(self, coordinates):
        self.coordinates = coordinates
    
    def setVelocities(self, velocities):
        self.velocities = velocities
    
    def getVelocities(self):
        return self.velocities
    
    def getCoordinateIds(self, boundaries):
        """
        returns the coordinateIDs of a segment with respect to 
        given boundaries
        """
        import lib.bin_classifier as bin_classifier
        return bin_classifier.getCoordinateIds(self.coordinates, boundaries)
        
    def __eq__(self, other_segment): 
        return self.__dict__ == other_segment.__dict__
    
    def isParent(self, segment):
        if segment.getParentIterationId() == self.getIterationId() and\
           segment.getParentBinId() == self.getBinId() and\
           segment.getParentSegmentId() == self.getId():
            return True
        return False

    
