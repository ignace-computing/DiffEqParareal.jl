This package provides an implementation of the parareal algorithm, which is a parallel-in-time method for the solution of time-dependent differential equations, such as ordinary differential equations (ODEs) and partial differential equations (PDEs). 

The algorithm requires two overlapping temporal grids, a coarse and a fine one. 
The algorithm operates iteratively: starting from a coarse numerical solution on the coarse grid, it iteratively refines this first guess, using updates that are computed *in parallel* in each coarse interval.

In its basic form, the parareal algorithm has the following form, where `u` is the solution of an ODE/PDE/...:

```math
u_{k+1}^{n+1} = C(u_{k+1}^n) +  F (u_k^n) - C(u_k^n) 
```

where `n` is the index of the coarse timepoints, `k` is the iteration number, and `C` and `F` are coarse and fine solution operators respectively.

TO DO: 
- [ ] Fix parallel execution. For now it only works correclty serially. 
- [ ] Once it works in parallel, test scalability.
- [ ] Add both interesting applications!
- [ ] Once all the above is done, add attractive applications, such as the Lorenz system.
