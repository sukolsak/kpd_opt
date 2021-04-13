mutable struct Input
    graph::DiGraph
    weights::Array{Float64,2}
    ndds::Array{Int64,1}
    max_cycle_length::Int64
    max_chain_length::Int64
    solver_instance::Any
    time_param::String
    time_factor::Float64
    start_time::Float64

    function Input(
        graph::DiGraph,
        weights::Array{Float64,2},
        ndds::Array{Int64,1},
        max_cycle_length::Int64,
        max_chain_length::Int64,
        solver_instance::Any,
        time_param::String,
        time_factor::Float64,
        start_time::Float64
    )
        new(graph, weights, ndds, max_cycle_length, max_chain_length, solver_instance, time_param, time_factor, start_time)
    end
end

mutable struct Output
    match_vertices::Array{Float64,1}
    match_edges::Any
    match_cycles::Union{Array{Float64,1}, Nothing}
    graph_cycles::Union{Array{Array{Int64,1},1}, Nothing}
    value::Float64
    optimal::Bool

    function Output()
        match_vertices = Float64[]
        match_edges = Float64[]
        match_cycles = Array{Float64}[]
        graph_cycles = Array{Float64}[]
        value = typemin(Float64)
        optimal = false
        new(match_vertices, match_edges, match_cycles, graph_cycles, value, optimal)
    end

    function Output(
        match_vertices::Array{Float64,1},
        match_edges::Any,
        match_cycles::Union{Array{Float64,1}, Nothing},
        graph_cycles::Union{Array{Array{Int64,1},1}, Nothing},
        value::Float64,
        optimal::Bool,
    )
        new(match_vertices, match_edges, match_cycles, graph_cycles, value, optimal)
    end
end
