using LightGraphsFlows
using SharedArrays
using JuMP
include("transport.jl")
const eps=1e-8

struct Intermediate
    m::Model
    x::JuMP.Containers.DenseAxisArray
    z::Array{VariableRef,1}
    in_flow::Array{VariableRef,1}
    out::Array{VariableRef,1}
    n::Int64
    max_cycle_length::Int64
    max_chain_length::Any
    graph::DiGraph
    graph_edges::Any
    index_to_edge::Dict{Any,Any}
    ndds::Array{Int64, 1}
    cb_graph::DiGraph
    caps::Array{Float64,2}
    graph_cycles::Any
    time_param::String
    time_factor::Float64
    start_time::Float64

    function Intermediate(input_data::Input)
        start_time = input_data.start_time
        graph, weights, ndds, max_cycle_length, max_chain_length, solver_instance = input_data.graph, input_data.weights, input_data.ndds, input_data.max_cycle_length, input_data.max_chain_length, input_data.solver_instance
        m = Model(solver_instance)
        n = nv(graph)
        index_to_edge = Dict()
        graph_edges = collect(edges(graph))
        for e in graph_edges
            index_to_edge[(e.src, e.dst)] = e
        end
        cycle_weights, c = get_cycles(graph, max_cycle_length, weights)
        num_cycles = size(c)[1]
        @variable(m, 1 >= x[graph_edges] >= 0, Int)
        @variable(m, 1 >= in_flow[1:n] >= 0, Int)
        @variable(m, 1 >= out[1:n] >= 0, Int)
        @variable(m, 1 >= z[1:num_cycles] >= 0, Int)

        @objective(m, Max, sum(x[Edge((i, j))] * weights[i, j] for i in 1:n, j in outneighbors(graph, i)) + sum(z[i] * cycle_weights[i] for i in 1:num_cycles))
        @constraint(m, c1[k=1:n], sum(x[Edge((i,k))] for i=inneighbors(graph, k)) == in_flow[k])
        @constraint(m, c2[k=1:n], sum(x[Edge((k,j))] for j=outneighbors(graph, k)) == out[k])
        @constraint(m, c3[k=setdiff(1:n, ndds)], out[k] <= in_flow[k])
        @constraint(m, c4[k=ndds], in_flow[k] == 0)
        @constraint(m, c5[k=setdiff(1:n, ndds)], in_flow[k] <= 1-sum(z[j] for j in cycles_part_of(c, k)))
        if length(ndds) == 0
            @constraint(m, d[e=graph_edges], x[e] == 0)
        end
        
        if max_chain_length < n-1
            @variable(m, 1 >= x_ndds[graph_edges, ndds] >= 0, Int)
            @variable(m, 1 >= in_ndds[1:n, ndds], Int)
            @variable(m, out_ndds[1:n, ndds], Int)

            @constraint(m, c6[e=graph_edges], sum(x_ndds[e, k] for k in ndds) == x[e])
            @constraint(m, c7[k=ndds], sum(x_ndds[e, k] for e in graph_edges) <= max_chain_length)
            @constraint(m, c8[k=1:n, ndd=ndds], sum(x_ndds[Edge((i,k)), ndd] for i=inneighbors(graph, k)) == in_ndds[k, ndd])
            @constraint(m, c9[k=1:n, ndd=ndds], sum(x_ndds[Edge((k,j)), ndd] for j=outneighbors(graph, k)) == out_ndds[k, ndd])
            @constraint(m, c10[k=setdiff(1:n, ndds), ndd=ndds], out_ndds[k, ndd] <= in_ndds[k, ndd])
        end
        cb_graph = deepcopy(graph)
        caps = zeros(Float64, n+1, n+1)
        add_vertices!(cb_graph, 1)
        for i in ndds
            add_edge!(cb_graph, n+1, i)
            caps[n+1, i] = 1
            caps[i, n+1] = Inf
        end
        println("Number of variables: ", length(all_variables(m)))
        # println("Number of constraints: ", num_constraints(m, Any, Any))
        time_param = input_data.time_param
        time_factor = input_data.time_factor
        new(m, x, z, in_flow, out, n, max_cycle_length, max_chain_length, graph, graph_edges, index_to_edge, ndds, cb_graph, caps, c, time_param, time_factor, start_time)
    end
end

function solve(problem::Intermediate, signal::SharedArray{Bool,1}=SharedArray{Bool}(1))
    output = Output()

    function general_callback(cb_data, MOI_instance)
        if length(problem.ndds) == 0
            return
        end
        if signal[1]
            println("Stopping solve")
            get_best_solution(problem, output, cb_data)
            con = @build_constraint(0 >= 1)
            MOI.submit(problem.m, MOI_instance(cb_data), con)
        end

        x_vals = callback_value.(Ref(cb_data), problem.x)
        in_vals = callback_value.(Ref(cb_data), problem.in_flow)
        pos_index = []
        for i in 1:problem.n
            if in_vals[i] > 0
                push!(pos_index, i)
            end
        end
        edges_sets = find_sets(x_vals, pos_index, problem)
        for pair in edges_sets
            i, edges_set = pair
            con = @build_constraint(sum(problem.x[e] for e in edges_set) >= problem.in_flow[i])
            MOI.submit(problem.m, MOI_instance(cb_data), con)
        end
    end

    function lazy_callback(cb_data)
        general_callback(cb_data, MOI.LazyConstraint)
    end
    MOI.set(problem.m, MOI.LazyConstraintCallback(), lazy_callback)

    function user_cut_callback(cb_data)
        general_callback(cb_data, MOI.UserCut)
    end
    MOI.set(problem.m, MOI.UserCutCallback(), user_cut_callback)

    tf = time()
    remaining_time = max(timeout-(tf-problem.start_time), 1)*problem.time_factor
    remaining_time = convert(Int32, round(remaining_time, digits=0))
    set_parameter(problem.m, problem.time_param, remaining_time)

    optimize!(problem.m)
    if termination_status(problem.m) == MOI.TIME_LIMIT
        println("Intermediate timed out")
        return nothing
    end
    return get_output(problem)
end

function get_output(problem::Intermediate)
    match_vertices = JuMP.value.(problem.in_flow)
    match_edges = JuMP.value.(problem.x)
    match_cycles = JuMP.value.(problem.z)
    graph_cycles = problem.graph_cycles
    value = JuMP.objective_value(problem.m)
    optimal = true
    output = Output(match_vertices, match_edges, match_cycles, graph_cycles, value, optimal)
    return output
end

function get_best_solution(problem::Intermediate, output::Output, cb_data::Any)
    return
end

"""
get_cycles
Returns list of weights for each cycle of length max k and list of said cycles.
"""
function get_cycles(graph::DiGraph, max_cycle_length::Int64, weights::Array{Float64,2})
    max_cycles = 1e8
    Ck = simplecycles_limited_length(graph, max_cycle_length, max_cycles)
    num_cycles = size(Ck)[1]
    println("Found $num_cycles cycles of length at most $max_cycle_length")
    w = zeros(Float64, num_cycles)
    for (i, cycle) in enumerate(Ck)
        n = size(cycle)[1]
        if n == 1
            w[i] = weights[cycle[1], cycle[1]]
        else
            w[i] = sum(weights[cycle[j], cycle[j+1]] for j in 1:n-1)
            w[i] += weights[cycle[n], cycle[1]]
        end
    end
    return (w, Ck)
end

"""
cycles_part_of
Returns list of indices corresponding to cycles that node k is a part of given
    the list C of cycles
"""
function cycles_part_of(C::Array{Array{Int64,1},1}, k::Int64)
    indices = zeros(Int64, 0)
    for (i,c) in enumerate(C)
        if k in c
            append!(indices, i)
        end
    end
    return indices
end

function find_sets(x_vals::JuMP.Containers.DenseAxisArray, pos_index::Array{Any, 1}, problem::Intermediate)
    n, ndds, graph, cb_graph, caps, index_to_edge = problem.n, problem.ndds, problem.graph, problem.cb_graph, problem.caps, problem.index_to_edge
    cb_caps = deepcopy(caps)
    for e in problem.graph_edges
        cb_caps[e.src,e.dst] = x_vals[e]
    end
    edges_sets = []
    for i in setdiff(pos_index, ndds)
        part1, part2, flow = minimum_cut(cb_graph, n+1, i, cb_caps)
        if flow <= 1-eps
            edges_set = []
            for j in setdiff(part2, ndds)
                vec = setdiff(inneighbors(graph, j),part2)
                for k in vec
                    push!(edges_set, index_to_edge[(k,j)])
                end
            end
            push!(edges_sets, (i, edges_set))
        end
    end
    return edges_sets
end

function minimum_cut(
        flow_graph::DiGraph,                 # the input graph
        source::Int64,                       # the source vertex
        target::Int64,                       # the target vertex
        capacity_matrix::Array{Float64,2},   # edge flow capacities
        # algorithm::AbstractFlowAlgorithm   # keyword argument for algorithm
    )
    flow, flow_matrix = maximum_flow(flow_graph, source, target, capacity_matrix, EdmondsKarpAlgorithm())
    # flow, flow_matrix = maximum_flow(flow_graph, source, target, capacity_matrix, PushRelabelAlgorithm())
    residual_matrix = zeros(nv(flow_graph),nv(flow_graph))
    for edge in edges(flow_graph)
        residual_matrix[edge.src,edge.dst] = max(0.0, capacity_matrix[edge.src,edge.dst] - flow_matrix[edge.src,edge.dst])
    end
    part1 = typeof(source)[]
    queue = [source]
    while !isempty(queue)
        node = pop!(queue)
        push!(part1, node)
        dests = [dst for dst in 1:nv(flow_graph) if residual_matrix[node,dst] > 0.0 && dst ∉ part1]
        append!(queue, dests)
    end
    part2 = [node for node in 1:nv(flow_graph) if node ∉ part1]
    return (part1, part2, flow)
end
