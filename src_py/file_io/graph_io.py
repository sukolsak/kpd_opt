import builtins
import itertools
import math
import os

import igraph as ig
import numpy as np

from match import infinity


def no_print(*args):
    return


def generate_graph(filespec, verbose=True):
    if not verbose:
        global print
        old_print = print
        print = no_print
    print()
    print("Creating graph from file")
    root_nodes, paired_nodes, terminal_nodes, edges, problem_data, forbidden_nodes = parse_input(filespec)
    ndds = list(range(len(root_nodes)))

    node_names = root_nodes + paired_nodes + terminal_nodes
    number_nodes = len(root_nodes) + len(paired_nodes) + len(terminal_nodes)

    id_map = {
        'r': 0,
        'p': len(root_nodes),
        't': len(root_nodes) + len(paired_nodes),
    }
    forbidden_nodes = [find_index(node, id_map) for node in forbidden_nodes]
    problem_data["forbiddenNodes"] = forbidden_nodes

    graph = ig.Graph(directed=True)
    graph.add_vertices(number_nodes)

    weights = np.zeros((number_nodes, number_nodes), dtype=float)
    ids_tuples = [None]*len(edges)
    for i, e in enumerate(edges):
        edge_id, origin, dest, edge_weight = e
        o_id = find_index(origin, id_map)
        d_id = find_index(dest, id_map)
        ids_tuples[i] = (o_id, d_id)
        weights[o_id, d_id] = float(edge_weight)
    graph.add_edges(ids_tuples)

    print()
    if not verbose:
        print = old_print
    return graph, weights, ndds, problem_data, node_names


def parse_input(file_spec):
    with open(file_spec) as f:
        lines = f.readlines()
        lines = list(map(str.strip, lines))
    problem_data_dict, i = read_problem_data_dict(lines)
    root_nodes, i = read_nodes(lines, i, "rootNodes")
    paired_nodes, i = read_nodes(lines, i, "pairedNodes")
    terminal_nodes, i = read_nodes(lines, i, "terminalNodes")
    try:
        forbidden_nodes, i = read_nodes(lines, i, "forbiddenNodes")
    except AssertionError as _:
        forbidden_nodes = []
    edges, i = read_edges(lines, i)
    problem_data_dict["rootNodes"] = len(root_nodes)
    problem_data_dict["pairedNodes"] = len(paired_nodes)
    problem_data_dict["terminalNodes"] = len(terminal_nodes)
    return root_nodes, paired_nodes, terminal_nodes, edges, problem_data_dict, forbidden_nodes


def read_nodes(nodes_list, line, str_format):
    nodes = []
    nodes_name, nodes_number = nodes_list[line].split(",")
    end_list = line + int(nodes_number)
    assert str.lower(nodes_list[end_list+1]) == str.lower("end" + str_format)
    print_property(nodes_name, nodes_number)
    for k in range(line, end_list):
        nodes += [nodes_list[k+1]]
    return nodes, end_list+2


def read_edges(edges_list, line):
    edges = []
    edges_name, edges_number = edges_list[line].split(",")
    print_property(edges_name, edges_number)
    end_list = line + int(edges_number)
    assert str.lower(edges_name) == "edges"
    for k in range(line+1, end_list+1):
        edges += [edges_list[k].split(",")]
    assert str.lower(edges_list[end_list+1]) == "endedges"
    return edges, end_list+2


def print_property(name, number):
    print("\t", name, number)


def find_index(node_name, id_map):
    letter = node_name[0]
    index = int(node_name[1:]) + id_map[letter] - 1
    return index


def read_problem_data_dict(str_list):
    problem_data_dict = dict()

    [chain_length_str, chain_length] = str_list[1].split(",")
    assert chain_length_str == "maxChainLength"
    problem_data_dict["max_chain_length"] = int(chain_length) if chain_length != "Infinity" else infinity

    [cycle_length_str, cycle_length] = str_list[2].split(",")
    assert cycle_length_str == "maxCycleLength"
    problem_data_dict["max_cycle_length"] = int(cycle_length) if cycle_length != "Infinity" else infinity

    cycle_bonus = str_list[3].split(",")
    assert cycle_bonus[0] == "cycleBonus"
    problem_data_dict["cycle_bonus"] = cycle_bonus[1]

    assert str_list[4] == "endProblemData"

    return problem_data_dict, 5


def sort_cycle(cycle, weights):
    n = len(cycle)
    min_vertex = math.inf
    min_index = 0
    for index, vertex in enumerate(cycle):
        if vertex < min_vertex:
            min_vertex = vertex
            min_index = index
    gen = [(i + min_index) % n for i in range(n)]
    cycle_sorted = [cycle[j] for j in gen]
    weights_sorted = [weights[j] for j in gen]
    return cycle_sorted, weights_sorted


def write_to_csv(cycle_chains_unsorted: list, filespec, node_names: list):
    cycle_chains_list = [None]*len(cycle_chains_unsorted)
    for i, cycle_chains in enumerate(cycle_chains_unsorted):
        order = cycle_chains[1]
        weights = cycle_chains[2]
        if cycle_chains[0] == 'cycle':
            order, weights = sort_cycle(order, weights)
        cycle_chains_list[i] = (cycle_chains[0], order, weights)
    cycle_chains_list.sort()

    print("Saving result to file:", os.path.realpath(filespec))
    mapping = {"chain": "endChain", "cycle": "endCycle"}
    file_str = ""
    for cycle_chains in cycle_chains_list:
        cat = cycle_chains[0]
        file_str += cat + "\n"
        list_vert = [node_names[i] for i in cycle_chains[1]]
        list_weights = cycle_chains[2]
        list_rows = [str(list_vert[i]) + "," + str(list_weights[i]) for i in range(len(list_vert))]
        str_list_vert = ", \n".join(list_rows)
        file_str += str_list_vert + "\n"
        file_str += mapping[cat] + "\n"
    with open(filespec, "w") as f:
        f.write(file_str)
