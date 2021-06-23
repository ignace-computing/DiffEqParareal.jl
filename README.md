This package provides an implementation of the parareal algorithm, which is a parallel-in-time method for the solution of time-dependent differential equations, such as ordinary differential equations (ODEs) and partial differential equations (PDEs). 

The algorithm requires two overlapping temporal grids, a coarse and a fine one. 
The algorithm operates iteratively: starting from a coarse numerical solution on the coarse grid, it iteratively refines this first guess, using updates that are computed *in parallel* in each coarse interval.

In its basic form, the parareal algorithm has the following form, where `u` is the solution of an ODE/PDE/...:

```math
u_{k+1}^{n+1} = C(u_{k+1}^n) +  F (u_k^n) - C(u_k^n) 
```

where `n` is the index of the coarse timepoints, `k` is the iteration number, and `C` and `F` are coarse and fine solution operators respectively.

## TO DO: 

- [ ] Make it run correctly in parallel. Best way forward is to run and read through`test_pararealParallel.jl`. In `pararealParallel.jl` I have added several suggestions for parallellism (using `Threads`, `Polyester` or `FLoops`, see lines 30-32), but none of them is working actually. I am not sure how we can resolve this.
- [ ] Test parallel scalability.
- [ ] Carry out experiments with interesting applications.
