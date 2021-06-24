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

MODEL = "simple"
if MODEL == "simple"
    function f(du,u,p,t)
        du[1] = -cos(u[1])*u[1]
    end
    u0 = [10.0]
    TSPAN = (0.0, 0.4)
    algC = Euler()
    algF = Euler()
    dtC = 0.1   # coarse timestep
    dtF = 0.01  # fine timestep
elseif MODEL == "Lorenz"
    function f(du,u,p,t)
        du[1] = 10.0(u[2]-u[1])
        du[2] = u[1]*(28.0-u[3]) - u[2]
        du[3] = u[1]*u[2] - (8/3)*u[3]
    end
    u0 = [1.0;0.0;0.0]
    TSPAN = (0.0, 100.0)
    algC = Euler()
    algF = Euler()
    dtC = 0.01   # coarse timestep
    dtF = 0.001  # fine timestep
end

println("Threads.nthreads() = ",Threads.nthreads())

# 1) reference solution
println("REFERENCE")
@btime begin
    probEXACT = ODEProblem(f, u0, TSPAN, dt=dtF)
    solEXACT = solve(probEXACT, algF, tstop=dtC)
end

probEXACT = ODEProblem(f, u0, TSPAN, dt=dtF)
solEXACT = solve(probEXACT, algF, tstop=dtC)

# 2) parareal solution
println("PARAREAL")
K = ceil(Int, (length(solEXACT.t)-1)/(Int(dtC/dtF)))

PP = Threads.nthreads()
@btime u_solution, sol = Parareal(f, u0, TSPAN, dtC, algC, dtF, algF, K=K, ncores=PP)
u_solution, sol = Parareal(f, u0, TSPAN, dtC, algC, dtF, algF, K=K, ncores=PP)

# @profilehtml u_solution, p2 = Parareal(f, u0, alg, TSPAN, dtC, dtF, ncores=PP)

# 3) plot the results
println("PLOT")
if length(u_solution) <= 101
    p2 = plot()
    plot!(p2, solEXACT.t, [z[1] for z in solEXACT.u], shape=:circle, markersize=2, label="exact")
    plot!(p2, sol.t, [z[1] for z in u_solution] , shape=:circle, markersize=2, label="parareal solution")
    plot!(p2, title="test parareal")
    plot!(p2, legend=:bottomright)
end

# 4) test the exactness property of parareal
println("TEST")
if length(u_solution) <= 11
    @testset "parareal exactness property" begin
        K = ceil(Int, (length(solEXACT.t)-1)/(Int(dtC/dtF)))
        for k=1:K+1
            u_solution, sol = Parareal(f, u0, TSPAN, dtC, algC, dtF, algF, K=K, ncores=PP)
            for n=1:k
                @test ≈(norm(solEXACT[1+(n-1)*Int(dtC/dtF)] - u_solution[n]), 0., atol=1e-12)
            end
        end
    end
end

println("END")
p2