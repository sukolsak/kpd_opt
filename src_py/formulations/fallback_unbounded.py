from lap import lapjv
import numpy as np

from formulations.formulation_abstract import Formulation
from utils.transport import Output, Input


class Unbounded(Formulation):
    def __init__(self, input_data: Input):
        self.graph = input_data.graph
        self.n = self.graph.vcount()
        self.can_run = (len(input_data.forbidden_nodes) == 0) \
                       and (input_data.cycle_length >= self.n) \
                       and ((input_data.chain_length >= self.n-1) or (len(input_data.ndds) == 0) or (input_data.chain_length == 0))
        self.has_chains = input_data.chain_length > 0
        self.input_data = input_data

    def solve(self):
        if not self.can_run:
            print("Unbounded fallback does not apply")
            return

        infty = 1e8
        n = self.n
        weight_matrix = np.full([n, n], infty)
        loops = set()
        for i in range(n):
            weight_matrix[i, i] = 0
        for e in self.graph.es():
            i = e.tuple[0]
            j = e.tuple[1]
            if i == j:
                loops.add(i)
            w = self.input_data.weights[i, j]
            weight_matrix[i, j] = -w
        ndds = self.input_data.ndds
        if self.has_chains:
            for j in ndds:
                for i in range(n):
                    if i in ndds:
                        continue
                    weight_matrix[i, j] = 0
        assignation_value, row_ind, col_ind = lapjv(weight_matrix)
        assignation_value = -assignation_value
        return self.get_output(assignation_value, row_ind, loops)

    def get_output(self, obj_val, row_ind, loops):
        match_edges = dict()
        value = 0
        for i, j in enumerate(row_ind):
            if j in self.input_data.ndds:  # ignore closure of cycles on ndds
                continue
            if (i == j) and i not in loops:  # ignore non assigned nodes
                continue
            match_edges[(i, j)] = 1
            value += self.input_data.weights[i, j]
        match_cycles = None
        graph_cycles = None
        output = Output(match_edges, match_cycles, graph_cycles, obj_val, True)
        return output
