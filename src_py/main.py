from file_io import graph_io
import sys

from match import chain_cycle_match, get_match_list, get_formulation, print_options

# TODO: check if not repeating cuts speeds up solves
# TODO: profile callbacks to make them faster
# TODO: move to graph-tool and implement cycle search in C++

from utils.transport import Input


def main(args):
    input_path = args[1]
    output_path = args[2]
    objective_fn = args[3]

    objective_fn = objective_fn.lower()

    print_options()
    print("Input path: %s" % input_path)
    print("Output path: %s" % output_path)
    print("Formulation: %s" % get_formulation(objective_fn).__name__)
    result = read_and_solve(input_path, objective_fn)
    if result is None:
        sys.exit(1)
    cycle_chains_list, formulation, input_data, output_data, node_names = result
    graph_io.write_to_csv(cycle_chains_list, output_path, node_names)


def read_and_solve(input_path, objective_fn, max_cycle_length=None, max_chain_length=None, decrease=True, verbose=True):
    graph, weights, ndds, problem_data, node_names = graph_io.generate_graph(input_path, verbose)

    if max_cycle_length is None:
        max_cycle_length = problem_data["max_cycle_length"]
    if max_chain_length is None:
        max_chain_length = problem_data["max_chain_length"]

    forbidden_nodes = problem_data["forbiddenNodes"]
    input_data = Input(graph, weights, ndds, max_cycle_length, max_chain_length, forbidden_nodes)
    result = chain_cycle_match(input_data, objective_fn, decrease)
    if result is None:
        return None
    output_data, formulation = result
    cycle_chains_list = get_match_list(output_data, graph, weights)
    return cycle_chains_list, formulation, input_data, output_data, node_names


if __name__ == '__main__':
    args = sys.argv
    main(args)
