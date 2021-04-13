using LightGraphs
include("match.jl")

"""
basic_MIP
Returns a JuMP model with the basic matching MIP. On its own, solves the problem with
unbounded chains and cycles.
"""
function parse_input(filespec::String)
    f = open(filespec)
    l = readlines(f)
    close(f)
    problem_data, i = read_problem_data(l)
    root_nodes, i = read_nodes(l, i, "rootNodes")
    paired_nodes, i = read_nodes(l, i, "pairedNodes")
    terminal_nodes, i = read_nodes(l, i, "terminalNodes")
    edges, i = read_edges(l, i)
    return root_nodes, paired_nodes, terminal_nodes, edges, problem_data
end

function read_nodes(list, line, format)
    nodes = []
    nodes_name, nodes_number = split(list[line], ",")
	num_elems = parse(Int64, nodes_number)
    end_list = line + num_elems
	if num_elems != 0
		for k in (line + 1):(end_list)
			push!(nodes, list[k])
		end
	else
		while list[end_list+1] == ""
			end_list += 1
		end
	end
   	@assert lowercase(list[end_list + 1]) == lowercase("end" * format)

    return nodes, end_list + 2
end

function read_edges(list, line)
    edges = []
    edges_name, edges_number = split(list[line], ",")
    end_list = line + parse(Int64, edges_number)
    @assert lowercase(edges_name) == "edges"
    for k in (line + 1):(end_list)
        push!(edges, (split(list[k], ",")))
    end
    @assert lowercase(list[end_list + 1]) == "endedges"
    return edges, end_list + 2
end

function generate_graph(filespec::String)
    input = parse_input(filespec)
    root_nodes, paired_nodes, terminal_nodes, edges, problem_data = input
    number_nodes = length(root_nodes) + length(paired_nodes) + length(terminal_nodes)
    ndds = 1:length(root_nodes)

    graph = LightGraphs.SimpleDiGraph(number_nodes)
    id_map = Dict('r' => 0,
                  'p' => length(root_nodes),
                  't' => length(root_nodes) + length(paired_nodes))
    weights = zeros(Float64, (number_nodes, number_nodes))
    for e in edges
        edge_id, origin, dest, edge_weight = e
        o_id = find_index(origin, id_map)
        d_id = find_index(dest, id_map)

        add_edge!(graph, o_id, d_id)
        weights[o_id, d_id] = parse(Float64, edge_weight)
    end
    return graph, weights, ndds, problem_data
end

function find_index(node_name, map)
    letter = node_name[1]
    index = parse(Int64, node_name[2:end]) + map[letter]
    return index
end
"""
Parses problem data from a text comma-separated file
Returns the problem data as a dictionnary and the line number to keep reading the file
"""
function read_problem_data(list)
    problem_data = Dict()

    chain_length_str, chain_length = split(list[2], ",")
    @assert chain_length_str == "maxChainLength"
    if chain_length == "Infinity"
        chain_length = typemax(Int64)
    else
        chain_length = parse(Int64, chain_length)
    end
    problem_data["max_chain_length"] = chain_length

    cycle_length_str, cycle_length = split(list[3], ",")
    @assert cycle_length_str == "maxCycleLength"
    if cycle_length == "Infinity"
        cycle_length = typemax(Int64)
    else
        cycle_length = parse(Int64, cycle_length)
    end
    problem_data["max_cycle_length"] = cycle_length

    cycle_bonus = split(list[4], ",")
    @assert cycle_bonus[1] == "cycleBonus"
    problem_data["cycle_bonus"] = cycle_bonus[2]

    @assert list[5] == "endProblemData"

    return problem_data, 6
end

function write_to_csv(cycle_chains_list, filespec)
    mapping = Dict("chain"=> "endChain", "cycle"=> "endCycle")
    str = ""
    for cycle_chains in cycle_chains_list
        cat = cycle_chains[1]
        str *= cat * "\n"
        list_vert = cycle_chains[2]
        str_list_vert = join(list_vert, ", \n")
        str *= str_list_vert * "\n"
        str *= mapping[cat] * "\n"
    end
    open(filespec, "w") do f
        write(f, str)
    end
end
