function Parareal(f, u0, alg, TSPAN, dtC, dtF; K=0, PARALLEL=true, ncores=1)

    # println("Parareal")

    # coarse time step
    prob1 = ODEProblem(f, u0, TSPAN, dt=dtC)
    alg1 = Euler()
    sol1 = solve(prob1, alg1, tstop=dtC)
    int1 = init(prob1, alg1, tstop=dtC, save_everystep=false)
    int3 = init(prob1, alg1, tstop=dtC, save_everystep=false)

    # fine time step
    prob2 = ODEProblem(f, u0, TSPAN, dt=dtF)
    alg2 = Euler()
    int2 = init(prob2, alg2, tstop=0.1, save_everystep=false)

    u_solution = [x[1] for x in sol1.u]
    u_parallel = [x[1] for x in sol1.u]
    u_coarse = [x[1] for x in sol1.u]

    if K == 0
        K = ceil(Int, (TSPAN[2]-TSPAN[1])/dtC)
    end
    for k=1:K

        # in the following loops, note that they start from k, exploiting the exactness property of parareal

        # PARALLEL LOOP
        TT= length(sol1.t)-1-k
        # @floop ThreadedEx(basesize=TT+ncores) for i = k:length(sol1.t)-1
        # @batch per=thread for i = k:length(sol1.t)-1
        # Threads.@threads for i = k:length(sol1.t)-1
        for i = k:length(sol1.t)-1
            DiffEqBase.set_ut!(int2, [copy(u_solution[i])], copy(sol1.t[i]))
            for _ in 1: ceil(Int, dtC/dtF)
                step!(int2)
            end
            u_parallel[i+1] = int2.u[1]
        end            

        # SEQUENTIAL LOOP
        for i = k:length(sol1.t)-1
            DiffEqBase.set_ut!(int1, [u_solution[i]], sol1.t[i])
            step!(int1)
            u_coarse[i+1] = int1.u[1]
            
            u_solution[i+1] = int1.u[1] + u_parallel[i+1] - u_coarse[i+1]
        end

        u_coarse = copy(u_solution)
    end

    return u_solution, sol1 
end
