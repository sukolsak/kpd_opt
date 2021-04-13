using LightGraphsMatching

struct Matching
    g::DiGraph
    w::Dict{Edge,Float64}
    solver_instance::Any
    function Matching(graph::DiGraph, weights::Array{Float64,2}, ndds::Array{Int64, 1}, solver_instance::Any)
        n = nv(graph)
        
        new(g, w, solver_instance)
    end
end

function solve(problem::Matching)
    match = maximum_weight_matching(g,with_optimizer(Cbc.Optimizer, logLevel=0),w)
    return problem.m, problem.x, problem.in_, nothing, nothing
end
