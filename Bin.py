#!/usr/bin/python3

class Bin(object):
    # points to the reference structure of this bin
    ref_coords = ""
    # the array of trajectories
    trajectories = []
    
    def __init__(self, ref_coords, trajectories):
        """ 
        @param ref_coords the path to reference coordinates defining the bin
        @param trajectories single or list of trajectories to 
               init the bin
        """
        self.ref_coords = ref_coords
        self.trajectories += trajectories
        
    def binProbability(self):
        """ Returns the cumulative probability of all bin trajectories"""
        pass
    
    def getNumberOfTrajectories(self):
        """Number of trajectories"""
        return len(trajectories)
    
    def forkTrajectories(self, N):
        """If the number of trajectories is smaller then N
           then we want to duplicated some trajectories 
           to fill the bin again"""
        pass
    
    def addTrajectory(self, trajectory):
        """
        Add the specified trajectory to this bin        
        @param """
        self.trajectories.append(trajectory)
    
    def delTrajectory(self, index):
        """
        Delete indexed trajectoy
        @param index of trajectory to delete
        """
        del self.trajectories[index]

