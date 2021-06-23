using DifferentialEquations
using LinearAlgebra, BenchmarkHistograms
using StatProfilerHTML
using Test
using Plots
plotly()

using Distributed
using FLoops
using Polyester


include("pararealParallel.jl")

function f(du,u,p,t)
    du[1] = -cos(u[1])*u[1]
end
u0 = [10.0]
TSPAN = (0.0, 0.4)



alg = Euler()
dtC = 0.1   # coarse timestep
dtF = 0.01  # fine timestep

println("Threads.nthreads() = ",Threads.nthreads())

# 1) reference solution
probEXACT = ODEProblem(f, u0, TSPAN, dt=dtF)
println("REFERENCE")
@btime solEXACT = solve(probEXACT, alg, tstop=dtC)
solEXACT = solve(probEXACT, alg, tstop=dtC)

# 2) parareal solution
K = Int((length(solEXACT.t)-1)//(Int(dtC/dtF)))
PP = Threads.nthreads()
println("PARAREAL")
@btime u_solution, sol = Parareal(f, u0, alg, TSPAN, dtC, dtF, K=K, ncores=PP)
u_solution, sol = Parareal(f, u0, alg, TSPAN, dtC, dtF, K=K, ncores=PP)

# @profilehtml u_solution, p2 = Parareal(f, u0, alg, TSPAN, dtC, dtF, ncores=PP)

# 3) plot the results
p2 = plot()
# plot!(p2, sol1.t, [z[1] fint1.u[1]or z in sol.u], shape=:circle, markersize=2, label="initial coarse")
plot!(p2, solEXACT.t, [z[1] for z in solEXACT.u], shape=:circle, markersize=2, label="exact")
if SEQUENTIAL
    plot!(p2, sol.t, u_solution , shape=:circle, markersize=2, label="parareal solution (sequential)")
end
plot!(p2, title="test Dahlquist "*string(alg))
plot!(p2, legend=:bottomright)

# 4) test the exactness property of parareal
# for now, the test only passes for the sequential version (PARALLEL=false)
println(" ")
@testset "test parareal exactness property" begin
    K = Int((length(solEXACT.t)-1)//(Int(dtC/dtF)))
    for k=1:K+1
        u_solution, sol = Parareal(f, u0, alg, TSPAN, dtC, dtF, K=K, ncores=PP)
        for n=1:k
            @test solEXACT[1+(n-1)*Int(dtC/dtF)][1] ≈ u_solution[n]
        end
    end
end

p2