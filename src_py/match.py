import sys
import time
import warnings

import igraph as ig
import numpy as np

from constants import timeout, solver_instance
from formulations.parallel import Parallel
from formulations.basic import Basic
from formulations.intermediate import Intermediate
from formulations.fallback import Fallback

from utils.graph_utils import find_recipients
from utils.transport import Output

infinity = int(sys.maxsize/2)
warnings.simplefilter('always', Warning)
formulations_dict = {"default": Parallel, "basic": Basic, "pctsp": Intermediate, "fallback": Fallback}


def get_formulation(objective_fn):
    return formulations_dict.get(objective_fn, formulations_dict["default"])


def print_options():
    print("Using solver:", solver_instance.get_solver_name())
    print("Timeout:", timeout, "seconds")


def get_str(length):
    return "Unbounded" if length >= infinity else str(length)


def chain_cycle_match(input_data, objective_fn, decrease=True, new_nodes=0, y=0):
    input_data.solver_instance = solver_instance
    chain_length = input_data.chain_length
    if input_data.cycle_length == 0:
        input_data.cycle_length = 1
    while True:
        cycle_length = input_data.cycle_length

        print('Maximum cycle length:', get_str(cycle_length))
        print('Maximum chain length:', get_str(input_data.chain_length))
        if new_nodes > 0:
            input_data.chain_length = chain_length+1
        if cycle_length == 0:
            time.sleep(0.5)
            message = "Invalid cycle length, aborting"
            warnings.warn(message, Warning, stacklevel=sys.maxsize)
            time.sleep(0.5)
            return None
        print()
        print("Initializing and solving formulation")
        t0 = time.perf_counter()
        input_data.start_time = t0
        formulation = get_formulation(objective_fn)
        problem = formulation(input_data)
        result = problem.solve()
        if result is None:
            if not decrease:
                return None

            if formulation == Fallback:
                return None

            input_data.cycle_length -= 1
            print()
            time.sleep(0.5)
            warnings.warn("Timed out, decreasing maximum cycle length to " + str(input_data.cycle_length), Warning, stacklevel=sys.maxsize)
            time.sleep(0.5)
            print()
        else:
            tf = time.perf_counter()
            if new_nodes > 0:
                postprocess_result(result[0], input_data.graph, new_nodes, y)
            print()
            print("Solution found")
            print("Solution time:", round(tf-t0, 3), "seconds")
            print()
            print("Objective value: %s" % round(result[0].value, 5))
            result[0].time = round(tf-t0, 7)
            break
    return result


def get_match_list(output_data: Output, graph: ig.Graph, weights: np.array):
    match_edges = output_data.match_edges
    match_cycles = output_data.match_cycles
    graph_cycles = output_data.graph_cycles

    n = graph.vcount()
    cycle_chains_list = []
    flag = [0]*n
    recipients = find_recipients(graph, match_edges)
    for base in range(n):
        if flag[base] == 0:
            current_vertex = base
            cycle_chain_length = 0
            cycle_chain = []
            cycle_chain_weights = []
            for k in range(n):
                next_vertex = recipients[current_vertex]
                if next_vertex == -1:  # no recipient, end of a chain
                    if cycle_chain_length > 0:
                        cycle_chain.append(current_vertex)
                        cycle_chain_weights.append(0)
                        cycle_chains_list.append(("chain", cycle_chain, cycle_chain_weights))
                    break
                if next_vertex == base:  # closed the cycle
                    cycle_chain.append(current_vertex)
                    cycle_chain_weights.append(weights[current_vertex, next_vertex])
                    cycle_chains_list.append(("cycle", cycle_chain, cycle_chain_weights))
                    break
                cycle_chain.append(current_vertex)
                cycle_chain_weights.append(weights[current_vertex, next_vertex])
                flag[next_vertex] = 1
                cycle_chain_length += 1
                current_vertex = next_vertex
    if match_cycles is not None:
        num_cycles = len(graph_cycles)
        for i in range(num_cycles):
            if match_cycles[i] > 1/2:
                cycle = graph_cycles[i]
                c = len(cycle)
                cycle_weights = [None]*c
                for j in range(c):
                    o = cycle[j]
                    d = cycle[(j+1) % c]
                    cycle_weights[j] = weights[o, d]
                cycle_chains_list.append(("cycle", cycle, cycle_weights))
    return cycle_chains_list
