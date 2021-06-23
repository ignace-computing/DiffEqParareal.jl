This package provides an implementation of the parareal algorithm, which is a parallel-in-time method for the solution of time-dependent differential equations, such as ordinary differential equations (ODEs) and partial differential equations (PDEs). 

The algorithm requires two overlapping temporal grids, a coarse and a fine one. 
The algorithm operates iteratively: starting from a coarse numerical solution on the coarse grid, it iteratively refines this first guess, using updates that are computed *in parallel* in each coarse interval, on a fine mesh.

In its basic form, the parareal algorithm has the following form, where $u$ is the solution of an ODE/PDE:

$$
u_{k+1}^{n+1} = \mathcal{C}_{\Delta t} (u_{k+1}^n) +  \mathcal{F}_{\Delta t} (u_k^n) - \mathcal{C}_{\Delta t}(u_k^n) 
$$

and where $n$ is the index of the coarse timepoints, $`k`$ is the iteration number, and $`C`$ and $`F`$ are coarse and fine solution operators to the differential equation that the parareal algorithm is aimed to solve (that depends of course on the specific application!).
