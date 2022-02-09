If you are interested in using this package, do not hesitate to contact me (see Github profile for contact information).

This package provides an implementation of the parareal algorithm, which is a parallel-in-time method for the solution of time-dependent differential equations, such as ordinary differential equations (ODEs) and partial differential equations (PDEs). Under the hood, it uses the [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl) package, that comes with a variety of solvers and possibilities. 

## The parareal algorithm.

The algorithm requires two overlapping temporal grids, a coarse and a fine one. 
The algorithm operates iteratively: starting from a coarse numerical solution on the coarse temporal grid, it iteratively refines this first guess, using updates that are computed *in parallel* in each coarse interval.

In its basic form, the parareal algorithm has the following form, where `u` is the solution of an ODE/PDE/...:

```math
u_{k+1}^{n+1} = C(u_{k+1}^n) +  F (u_k^n) - C(u_k^n) 
```

where `n` is the index of the coarse timepoints, `k` is the iteration number, and `C` and `F` are coarse and fine solution operators, respectively. This algorithm was first proposed in [1]

## TO DO: 

- [ ] Improve parallel scalability! For now scaling results are awful, parallel execution takes longer than serial execution.
- [ ] Carry out experiments with interesting applications.
- [x] Make the implemented algorithm run correctly in parallel. Best way forward is to run and read through`test_pararealParallel.jl`. In `pararealParallel.jl` I have added several suggestions for parallellism (using [`Threads.jl`](https://docs.julialang.org/en/v1/base/multi-threading/), [`Polyester.jl`](https://github.com/JuliaSIMD/Polyester.jl) or [`FLoops.jl`](https://github.com/JuliaFolds/FLoops.jl), see lines 30-32), but none of them is working actually. I am not sure how we can resolve this.

## References
[1] J.-L. Lions, Y. Maday, and G. Turinici, “A ‘parareal’ in time discretization of PDE’s,” Comptes Rendus de l’Académie des Sciences - Series I - Mathematics, vol. 332, pp. 661–668, 2001 [Online]. Available at: http://dx.doi.org/10.1016/S0764-4442(00)01793-6
