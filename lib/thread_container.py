#!/usr/bin/env python2
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

class ThreadContainer:
    """
    A simple container class that manages bunches of threads.
    """
    def __init__(self):
        self.jobs = []

    def runJobs(self):
        """
        Run a list of jobs in parallel
        @param job_list of jobs to run
        @return a new empty job_list
        """
        # Start all jobs in parallel
        for job in self.jobs:
            job.start()
        # Wait until threads are finished
        for job in self.jobs:
            job.join()
        # Reset the job list
        self.jobs = []
    
    def appendJob(self, thread):
        """
        Append a job to the queue
        """
        self.jobs.append(thread)

    def getNumberOfJobs(self):
        """
        @return how many jobs are in the list
        """
        return len(self.jobs)
    
    def emptyQueue(self):
        """
        Remove all jobs from the queue without
        running them.
        """
        self.jobs = []
