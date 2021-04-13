using Distributed

dt = 0.1
exit_time = 1

n_formulations = 2

max_iterations = max_cycle_length
iteration = 0
results = nothing
finished = nothing
wids = nothing

function clean_results()
    global results
    global finished
    results = Vector{Any}(nothing, n_formulations)
    finished = Array{Bool,1}(undef, n_formulations)
end

function create_procs()
    extra_procs = max(max_iterations*n_formulations+1-nprocs(), 0)
    addprocs(extra_procs)
    proc_min = 2+iteration*n_formulations
    proc_max = proc_min+n_formulations-1
    global wids
    wids = collect(proc_min:proc_max)
    @everywhere include("includes.jl")
end

clean_results()
create_procs()

struct Parallel
    args::Any

    function Parallel(args...)
        new(args)
    end
end

@everywhere function solve_basic(index, args)
    println("Basic started")
    instance = Basic(args...)
    result = solve(instance)
    println("Basic finished")
    @spawnat 1 begin
        results[index] = result
        finished[index] = true
    end
end

@everywhere function solve_intermediate(index, args)
    println("Intermediate started")
    instance = Intermediate(args...)
    result = solve(instance)
    println("Intermediate finished")
    @spawnat 1 begin
        results[index] = result
        finished[index] = true
    end
end

function get_first(args)
    implementations = [solve_basic, solve_intermediate]

    clean_results()
    create_procs()
    n = n_formulations
    for i in 1:n
        wid = wids[i]
        f = Future(wid)
        @spawnat wid implementations[i](i, args)
    end
    first_finished = false
    result = nothing
    while !first_finished
        first_finished = all_finished()
        for i in 1:n
            if !isnothing(results[i])
                result = results[i]
                first_finished = true
                break
            end
        end
        sleep(dt)
    end
    exit_processes()
    global iteration
    iteration += 1
    return result
end

function all_finished()
    return all(finished)
end

function exit_processes()
    for i in wids
        @async interrupt(i)
    end
    sleep(exit_time)
end

function solve(problem::Parallel)
    return get_first(problem.args)
end
