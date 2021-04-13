from formulations.formulation_abstract import Formulation
from utils.max_weight_matching import max_weight_matching
from utils.transport import Input, Output

import igraph as ig


class Matching(Formulation):
    def __init__(self, input_data: Input):
        self.can_run = (len(input_data.forbidden_nodes) == 0) \
                       and (input_data.cycle_length == 2) \
                       and (len(input_data.ndds) == 0 or input_data.chain_length <= 1)
        self.input_data = input_data

    def solve(self):
        if not self.can_run:
            print("Matching fallback does not apply")
            return
        graph_orig = self.input_data.graph
        graph = graph_orig.as_undirected(mode="mutual")
        if self.input_data.chain_length == 1:
            ndd_edges = []
            for n in self.input_data.ndds:
                for v in graph_orig.neighbors(n, ig.OUT):
                    ndd_edges.append((n, v))
            graph.add_edges(ndd_edges)

        weights_orig = self.input_data.weights
        weights = dict()
        n = graph.vcount()
        m = graph.ecount()
        repeated_vertices = {}

        for edge in graph.es:
            pair = edge.tuple
            w1 = weights_orig[pair[0], pair[1]]
            w2 = weights_orig[pair[1], pair[0]]
            w = w1 + w2
            weights[pair] = w
            edge["weight"] = w
            if pair[0] == pair[1]:
                w /= 2
                repeated_vertices[pair[0]] = w

        repeated_list = list(repeated_vertices.keys())
        repeated_list.sort()
        new_edges = [(v, n+i) for i, v in enumerate(repeated_list)]
        new_weights = [repeated_vertices[v] for v in repeated_list]
        graph.add_vertices(len(repeated_list))
        graph.add_edges(new_edges)
        graph.es[m:]["weight"] = new_weights

        for pair in zip(new_edges, new_weights):
            v = pair[0]
            w = pair[1]
            weights[v] = w

        result = max_weight_matching(graph)
        return self.get_output(result, weights, n)

    @staticmethod
    def get_output(result, weights, n):
        match_edges = dict()
        obj_val = 0

        for key, value in result.items():
            pair = (key, value)
            pair_rev = tuple(reversed(pair))
            w1 = weights.get(pair, 0)
            w2 = weights.get(pair_rev, 0)
            w = w1 + w2
            obj_val += w
            if key >= n or value >= n:
                v = min(key, value)
                pair = (v, v)
            match_edges[pair] = 1
        obj_val /= 2

        match_cycles = None
        graph_cycles = None
        output = Output(match_edges, match_cycles, graph_cycles, obj_val, True)
        return output
