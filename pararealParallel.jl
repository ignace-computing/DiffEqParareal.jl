function Parareal(f, u0, TSPAN, dtC, algC, dtF, algF; K=0, ncores=1)

    # coarse time step
    prob1 = ODEProblem(f, u0, TSPAN, dt=dtC)
    sol1 = solve(prob1, algC, tstop=dtC)
    int1 = init(prob1, algC, tstop=dtC, save_everystep=false)
    N = length(sol1.t)-1

    # fine time step
    prob_fine = ODEProblem(f, u0, TSPAN, dt=dtF)
    parallel_ints = [init(prob_fine, algF, tstop=0.1, save_everystep=false) for _ in 1:ncores]

    u_solution = [x for x in sol1.u]
    u_parallel = [x for x in sol1.u]
    u_coarse = [x for x in sol1.u]

    if K == 0
        K = ceil(Int, (TSPAN[2]-TSPAN[1])/dtC)
    end
    for k=1:K

        # in the following loops, note that they start from k, exploiting the exactness property of parareal

        # PARALLEL LOOP
        TT= length(sol1.t)-1-k
        # @floop ThreadedEx(basesize=TT+ncores) for i = k:length(sol1.t)-1
        # @batch per=thread for i = k:length(sol1.t)-1
        Threads.@threads for i = k:length(sol1.t)-1
        # for i = k:length(sol1.t)-1
            int_p = parallel_ints[Threads.threadid()]
            DiffEqBase.set_ut!(int_p, copy(u_solution[i]), sol1.t[i])
            for _ in 1: ceil(Int, dtC/dtF)
                step!(int_p)
            end
            u_parallel[i+1] = int_p.u
        end            

        # SEQUENTIAL LOOP
        for i = k:length(sol1.t)-1
            DiffEqBase.set_ut!(int1, copy(u_solution[i]), sol1.t[i])
            step!(int1)
            u_coarse[i+1] = int1.u
            
            u_solution[i+1] = int1.u + u_parallel[i+1] - u_coarse[i+1]
        end

        u_coarse = copy(u_solution)
    end

    return u_solution, sol1 
end
