# hdWE
**hdWE** is an implementation of the Weighted Ensemble methodology originally
invented by Huber and Kim in 1996.
**hdWE** stands for **Hyper-Dimensional Weighted Ensemble**.

The code is purely written in Python 2 to guarantee a maximum of platform
independence and provide an easy start for coding newcomers who want to understand
and expand the implementation.
The WE algorithm is subdivided in two basic parts: While the bookkeeping of bin
structures, WE iterations, analysis of bin coordinates, and  stochastic weights
of the trajectories is handled by the hdWE program itself, the molecular dynamics
simulations are outsourced to an external MD software suite (e.g. AMBER, Gromacs).
This subdivision allowed us to focus our efforts on the design and performance
improvement of the WE implementation itself while taking advantage of state
of the art MD implementations at the same time.
Two key benefits of this approach should be noted.
First, as the weighted ensemble does not require modifications of the Hamiltonian,
we rely on the fastest available implementations of pure MD algorithms which are
nowadays typically accelerated with the help of General Purpose graphics processing
Units (GPUs) and allow us to access timescale of several microseconds for a typical
biomolecular system of approximately 50000 atoms.
Second, during the propagation step of a WE iteration, the trajectories are
completely independent from each other and can therefore be propagated in
multiple parallel MD threads for the next WE time step.
This form of trivial parallelizability is an intrinsic advantage of the WE algorithm.
After the trajectories have been propagated, their end structures are evaluated
with respect to their new position on the binning coordinates and resorted into the bins.
Then the trajectories are split and merged according to the resampling mechanism
to ensure that the bins are evenly filled.
Eventually the algorithm reenters the MD propagation step.
