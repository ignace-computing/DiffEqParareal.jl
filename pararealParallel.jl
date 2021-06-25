function Parareal(f, u0, TSPAN, dtC, algC, dtF, algF; K=0, ncores=1)
  # coarse time step
  prob1 = ODEProblem(f, u0, TSPAN, dt=dtC);
  prob_fine = ODEProblem(f, u0, TSPAN, dt=dtF);
  Parareal(prob1, prob_fine, u0, TSPAN, dtC, algC, dtF, algF, K, ncores)
end
ceilint(x) = ceil(Int, x)
function Parareal(prob1, prob_fine, u0, TSPAN, dtC, algC, dtF, algF, K, ncores)
  # t = time_ns()
  # parpart = zero(t)
    sol1 = solve(prob1, algC, tstop=dtC)
    int1 = init(prob1, algC, tstop=dtC, save_everystep=false)
    N = length(sol1.t)-1

    # fine time step
    parallel_ints = [init(prob_fine, algF, tstop=0.1, save_everystep=false) for _ in 1:ncores]

    u_solution = [x for x in sol1.u]
    u_parallel = [x for x in sol1.u]
    u_coarse = [x for x in sol1.u]
    u_solution_buffers = similar(u_solution, (Sys.CPU_THREADS)::Int)
    let x = first(u_solution)
      for i ∈ eachindex(u_solution_buffers)
        u_solution_buffers[i] = similar(x)
      end
    end
    u_sol_buf0 = u_solution_buffers[1]
    if K == 0
        K = ceil(Int, (TSPAN[2]-TSPAN[1])/dtC)
    end
    for k=1:K

        # in the following loops, note that they start from k, exploiting the exactness property of parareal

        # PARALLEL LOOP
      TT= length(sol1.t)-1-k
      # part_start = time_ns()
        #@floop ThreadedEx(basesize=TT+ncores) for i = k:length(sol1.t)-1
        # @batch per=core for i = k:length(sol1.t)-1
        #Threads.@threads for i = k:length(sol1.t)-1
        for i = k:length(sol1.t)-1
            tid = Threads.threadid()
            u_sol_buff = copyto!(u_solution_buffers[tid], u_solution[i])
            int_p = parallel_ints[tid]
            DiffEqBase.set_ut!(int_p, u_sol_buff, sol1.t[i])
            for _ in 1: ceilint(dtC/dtF)
                step!(int_p)
            end
            u_parallel[i+1] .= int_p.u
        end
      # parpart += time_ns() - part_start

        # SEQUENTIAL LOOP
        for i = k:length(sol1.t)-1
            DiffEqBase.set_ut!(int1, copyto!(u_sol_buf0, u_solution[i]), sol1.t[i])
            step!(int1)
            u_coarse[i+1] = int1.u
            
            @. u_solution[i+1] = int1.u + u_parallel[i+1] - u_coarse[i+1]
        end

        copyto!(u_coarse, u_solution)
        # u_coarse = u_solution
    end

    # return u_solution, sol1, parpart / (time_ns() - t)
    return u_solution, sol1
end
