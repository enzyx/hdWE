"""
A helper module to pack the cleanup process into threads and 
do not block the computationally way more intensive simulation
flow. 
"""
import threading

class Cleanup(object):
    """
    Contains a list of cleanup processes and stores
    data connected to cleanup.
    """
    def __init__(self, md_module, number_of_threads, compress_iteration, 
                 compress_closest_mask, keep_coords_frequency, keep_coords_segments, debug=False):
        self.cleanup_threads       = []
        self.md_module             = md_module
        self.compress_iteration    = compress_iteration
        self.compress_closest_mask = compress_closest_mask
        self.keep_coords_frequency = keep_coords_frequency
        self.keep_coords_segments  = keep_coords_segments
        self.debug                 = debug
        self.number_of_threads     = number_of_threads
        self.is_threaded           = bool(number_of_threads > 1)

    def _cleanup(self, iteration):
        """
        Internal cleanup procedure. Should not be called from 
        outside this module. use self.doCleanup() instead
        """
        if self.compress_iteration:
            self.md_module.compressIteration(iteration, self.compress_closest_mask)
        if  self.keep_coords_frequency == 0 or iteration.getId() % self.keep_coords_frequency != 0:
            self.md_module.removeCoordinateFiles(iteration, self.compress_iteration, self.keep_coords_segments)
        if self.debug:
            print("Cleanup process done (iteration {})".format(iteration.getNameString()))

    def doCleanup(self, iteration):
        if self.is_threaded:
            cleanup_thread = threading.Thread(target=self._cleanup, args=(iteration,))    
            cleanup_thread.start()
            self.cleanup_threads.append(cleanup_thread)
        else:
            self._cleanup(iteration)
        
        # Avoid to accumulate cleanup processes by waiting for 
        # them to finish every N iteration but do not wait for 
        # the very last one which was just submitted.
        if self.is_threaded and (len(self.cleanup_threads) >= self.number_of_threads):
            if self.debug:
                print("Waiting for cleanup threads to finish")
            for cleanup_thread in self.cleanup_threads[:-1]:
                cleanup_thread.join()
            self.cleanup_threads = self.cleanup_threads[-1:]
    
    def finishCleanupThreads(self):
        for cleanup_thread in self.cleanup_threads:
            cleanup_thread.join()
    