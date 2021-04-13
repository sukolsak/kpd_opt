using ArgParse

include("constants.jl")
include("match.jl")
include("io.jl")
include("transport.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "input_path"
            help = "Path to the input Data containing the graph to match."
            required = true
        "output_path"
            help = "Path where the output will be writen."
            required = true
        "objective_fn"
            help = "Type of objective function (basic, pctsp or parallel)"
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Input path: $(parsed_args["input_path"])")
    println("Output path: $(parsed_args["output_path"])")
    println("Objective function: $(parsed_args["objective_fn"])")
    println("Timeout: $timeout seconds")
    println()
    graph, weights, ndds, problem_data = generate_graph(parsed_args["input_path"])
    max_cycle_length = problem_data["max_cycle_length"]
    max_chain_length = problem_data["max_chain_length"]
    
    # Override lengths:
    # max_cycle_length = 2
    # max_chain_length = 1
    # max_cycle_length = typemax(Int64)
    # max_chain_length = typemax(Int64)
    
    input_data = Input(graph, weights, Array(ndds), max_cycle_length, max_chain_length, nothing, time_param, time_factor, 0.0)
    match_vertices, match_edges, match_cycles, graph_cycles, value = chain_cycle_match(input_data, parsed_args["objective_fn"])
    cycle_chains_list = get_match_list(match_edges, match_cycles, graph_cycles, graph)
    write_to_csv(cycle_chains_list, parsed_args["output_path"])
end

main()
