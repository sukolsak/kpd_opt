using SharedArrays
include("transport.jl")
const eps=1e-8

struct Basic
    m::Model
    x::JuMP.Containers.DenseAxisArray
    z::JuMP.Containers.DenseAxisArray
    in_flow::Array{VariableRef,1}
    out_flow::Array{VariableRef,1}
    n::Int64
    max_cycle_length::Int64
    max_chain_length::Int64
    graph::DiGraph
    time_param::String
    time_factor::Float64
    start_time::Float64

    function Basic(input_data::Input)
        # start_time = time()
        start_time = input_data.start_time
        graph, weights, ndds, max_cycle_length, max_chain_length, solver_instance = input_data.graph, input_data.weights, input_data.ndds, input_data.max_cycle_length, input_data.max_chain_length, input_data.solver_instance
        m = Model(solver_instance)
        n = nv(graph)
        graph_edges = collect(edges(graph))
        patients = setdiff(1:n, ndds)

        @variable(m, 1 >= x[graph_edges] >= 0, Int)
        @variable(m, 1 >= y[graph_edges] >= 0, Int)
        @variable(m, z[graph_edges] >= 0, Int)
        @variable(m, 1 >= in_flow[1:n] >= 0, Int)
        @variable(m, 1 >= out_flow[1:n] >= 0, Int)

        @objective(m, Max, sum(z[e]*weights[e.src, e.dst] for e in graph_edges))
        @constraint(m, c0[e=graph_edges], x[e]+y[e] == z[e])
        @constraint(m, c1[k=1:n], sum(x[Edge((i,k))] for i=inneighbors(graph, k)) == in_flow[k])
        @constraint(m, c2[k=1:n], sum(x[Edge((k,j))] for j=outneighbors(graph, k)) == out_flow[k])
        @constraint(m, c3[k=patients], out_flow[k] <= in_flow[k])
        @constraint(m, c4[k=ndds], in_flow[k] == 0)

        if max_chain_length < n-1
            @constraint(m, c5[p=patients], out_flow[p] == in_flow[p])
            @constraint(m, c6[k=ndds], out_flow[k] == 0)
        
            @variable(m, 1>= x_ndds[ndds, graph_edges] >= 0, Int)
            @variable(m, 1>= in_ndds[ndds, 1:n] >= 0, Int)
            @variable(m, 1>= out_ndds[ndds, 1:n] >= 0, Int)
        
            @constraint(m, c7[k=ndds], sum(x_ndds[k, e] for e in graph_edges) <= max_chain_length)
            @constraint(m, c8[v=ndds, k=1:n], sum(x_ndds[v, Edge((i,k))] for i=inneighbors(graph, k)) == in_ndds[v, k])
            @constraint(m, c9[v=ndds, k=1:n], sum(x_ndds[v, Edge((k,j))] for j=outneighbors(graph, k)) == out_ndds[v, k])
        
            @constraint(m, c10[v=ndds, p=patients], out_ndds[v, p] <= in_ndds[v, p])
            @constraint(m, c11[k=1:n], in_flow[k] + sum(in_ndds[v, k] for v in ndds) <= 1)
            @constraint(m, c13[k=1:n], out_flow[k] + sum(out_ndds[v, k] for v in ndds) <= 1)
            @constraint(m, c14[v=ndds, k=ndds], in_ndds[v, k] == 0)
        
            @constraint(m, c12[e=graph_edges], y[e] == sum(x_ndds[v, e] for v in ndds))
        else
            @constraint(m, c12[e=graph_edges], y[e] == 0)
        end

        time_param = input_data.time_param
        time_factor = input_data.time_factor
        new(m, x, z, in_flow, out_flow, n, max_cycle_length, max_chain_length, graph, time_param, time_factor, start_time)
    end
end

function solve(problem::Basic, signal::SharedArray{Bool,1}=SharedArray{Bool}(1))
    found = true
    cycles = []
    output = Output()
    while found
        if signal[1]
            println("Stopping solve")
            get_best_solution(problem, output)
            return output
        end
        for element in cycles
            cycle_length = element[1]
            cycle = element[2]
            # println(cycle_length, " - ", cycle)
            # @constraint(problem.m, sum(problem.z[e] for e in cycle) <= cycle_length-1)
            @constraint(problem.m, sum(problem.z[e] for e in cycle) <= cycle_length-1)
        end
        tf = time()
        remaining_time = max(timeout-(tf-problem.start_time), 1)*problem.time_factor
        remaining_time = convert(Int32, round(remaining_time, digits=0))
        set_parameter(problem.m, problem.time_param, remaining_time)
        optimize!(problem.m)
        # println(JuMP.objective_value(problem.m))
        if termination_status(problem.m) == MOI.TIME_LIMIT
            println("Basic timed out")
            return nothing
        end
        z_val = JuMP.value.(problem.z)
        found, cycles = check_cycle_lengths(z_val, problem.n, problem.max_cycle_length, problem.max_chain_length, problem.graph, true)  # only find the first violating cycle, seems to be faster
        # found, cycles = check_cycle_lengths(z_val, problem.n, problem.max_cycle_length, problem.max_chain_length, problem.graph, false)  # find all violating cycles in the current solution
    end
    return get_output(problem)
end

function get_output(problem::Basic)
    match_vertices = JuMP.value.(problem.in_flow)  # todo: correct to total in_flow
    match_edges = JuMP.value.(problem.z)
    match_cycles = nothing
    graph_cycles = nothing
    value = JuMP.objective_value(problem.m)
    optimal = true
    output = Output(match_vertices, match_edges, match_cycles, graph_cycles, value, optimal)
    return output
end

function get_best_solution(problem::Basic, output::Output)
    return
end

function print_values(x_val)
    for key in keys(x_val)
        value = x_val[key]
        if value > eps
            println(key, value)
        end
    end
end

function check_cycle_lengths(x_val, n::Int64, max_cycle_length::Int64, max_chain_length::Int64, graph::DiGraph, early_stop::Bool=false)
    # print_values(x_val)
    flag = zeros(Int64, n)
    cycles = []
    for i in 1:n
        if flag[i] != 0
            continue
        end
        current_vertex = i
        cycle_length = 0
        cycle = []
        while true
            flag[current_vertex] = 1
            next_vertex = find_recipient(current_vertex, x_val, graph)
            if next_vertex == 0 # no recipient, end of a chain
                if cycle_length > max_chain_length
                    push!(cycles, (cycle_length, cycle))
                    if early_stop
                        return true, cycles
                    end
                end
                break
            end
            push!(cycle, Edge((current_vertex, next_vertex)))
            flag[next_vertex] = 1
            cycle_length += 1
            if next_vertex == i # found a cycle
                if cycle_length > max_cycle_length
                    push!(cycles, (cycle_length, cycle))
                    if early_stop
                        return true, cycles
                    end
                end
                break
            end
            current_vertex = next_vertex
        end
    end
    return length(cycles)>0, cycles
end

function find_recipient(i::Int64, x_val, graph::DiGraph)
    for j in outneighbors(graph, i)
        if x_val[Edge((i,j))] > 0.5
            return j
        end
    end
    return 0
end
