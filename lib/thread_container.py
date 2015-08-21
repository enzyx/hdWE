#!/usr/bin/env python2

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
