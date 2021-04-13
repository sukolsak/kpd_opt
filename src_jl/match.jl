using JuMP, LightGraphs, GLPK, GLPKMathProgInterface

time_factor = 1.0

solver_instance = with_optimizer(GLPK.Optimizer, tm_lim=timeout*1000); time_param="tm_lim"; time_factor=1000.0
# using Gurobi; const GRB_ENV = Gurobi.Env(); solver_instance = with_optimizer(Gurobi.Optimizer, GRB_ENV, OutputFlag=0, LazyConstraints=1, TimeLimit=timeout); time_param = "TimeLimit"
# using CPLEX, MathProgBase; solver_instance = with_optimizer(CPLEX.Optimizer, CPX_PARAM_TILIM=timeout)

include("basic.jl")
include("intermediate.jl")
include("parallel.jl")
include("transport.jl")

"""
chain_cycle_match
Returns a maximum-weight matching in the form of three elements:
- An array of matched vertices
- An array of matched edges
- The corresponding match value.
"""
function chain_cycle_match(input_data::Input, objective_fn::String = "basic")
    if input_data.max_cycle_length == 0
        input_data.max_cycle_length = 1
    end
    while true
        t0 = time()
        @time begin
            if input_data.max_cycle_length <= 0 || input_data.max_chain_length < 0
                println("Error in parameters:")
                println("\tMaximum cycle length: $(input_data.max_cycle_length)")
                println("\tMaximum chain length: $(input_data.max_chain_length)")
                exit(1)
            end
            input_data.solver_instance = solver_instance
            input_data.time_param = time_param
            input_data.time_factor = time_factor
            input_data.start_time = t0

            println("Solving for parameters:")
            println("\tMaximum cycle length: $(input_data.max_cycle_length)")
            println("\tMaximum chain length: $(input_data.max_chain_length)")
            
            if objective_fn == "basic"
                problem = Basic(input_data)
            elseif objective_fn == "pctsp"
                problem = Intermediate(input_data)
            else
                problem = Parallel(input_data)
            end

            output_data = solve(problem)
        end

        if !isnothing(output_data) && output_data.optimal
            match_vertices = output_data.match_vertices
            match_edges = output_data.match_edges
            match_cycles = output_data.match_cycles
            graph_cycles = output_data.graph_cycles
            value = round(output_data.value, digits=5)
            
            println("\nObjective value: $value")
            return (match_vertices, match_edges, match_cycles, graph_cycles, value)
        end
        println()
        println("Solution is not optimal, decreasing maximum lengths")
        input_data.max_cycle_length = input_data.max_cycle_length-1
    end
end

"""
generate_graph
Creates a LightGraph directed graph from an adjacency edge-weight matrix.
"""
function generate_graph(weights::Array{Float64,2})
    @assert size(weights)[1] == size(weights)[2]
    n = size(weights)[1]
    graph = DiGraph(n)
    for i in 1:n
        for j in 1:n
            if weights[i,j] > 0
                add_edge!(graph, i, j)
            end
        end
    end
    return graph
end

function get_match_list(match_edges, match_cycles, graph_cycles, graph::DiGraph)
    n = nv(graph)
    cycle_chains_list = []
    flag = zeros(Int64, n)
    for i in 1:n
        if flag[i] == 0
            current_vertex = i
            cycle_chain_length = 0
            cycle_chain = []
            for k in 1:n  # replaces while true in case the match_edges is not a valid match.
                next_vertex = find_recipient(current_vertex, match_edges, graph)
                if next_vertex == 0 # no recipient, end of a chain
                    if cycle_chain_length > 0
                        push!(cycle_chain, current_vertex)
                        push!(cycle_chains_list, ("chain", cycle_chain))
                    end
                    break
                end

                if next_vertex == i  # closed the cycle
                    push!(cycle_chain, current_vertex)
                    push!(cycle_chains_list, ("cycle", cycle_chain))
                    break
                end
                push!(cycle_chain, current_vertex)
                flag[next_vertex] = 1
                cycle_chain_length += 1
                current_vertex = next_vertex
            end
        end
    end
    if !isnothing(match_cycles)
        num_cycles = size(graph_cycles)[1]
        for i in 1:num_cycles
            if match_cycles[i] > 0.5
                push!(cycle_chains_list, ("cycle", graph_cycles[i]))
            end
        end
    end
    return cycle_chains_list
end
