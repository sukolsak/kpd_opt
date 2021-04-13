import math

import igraph as ig

from constants import eps
from formulations.formulation_abstract import Formulation
from utils.graph_utils import simple_cycles_limited_length, simple_cycles_parallel
from utils.transport import Input, Output
from utils.utils import get_remaining_time, setdiff


class Intermediate(Formulation):
    def __init__(self, input_data: Input):
        self.start_time = input_data.start_time
        self.timed_out = False

        self.m = input_data.solver_instance("basic")
        self.max_cycle_length = input_data.cycle_length
        self.max_chain_length = input_data.chain_length
        self.graph = graph = input_data.graph

        self.ndds = input_data.ndds
        self.weights = input_data.weights

        self.forbidden_nodes = input_data.forbidden_nodes

        n = graph.vcount()
        self.n = n
        vertices = graph.vs
        edges = graph.es
        self.edge_tuples = [e.tuple for e in edges]
        self.vertex_indices = [v.index for v in vertices]
        self.set_neighbors()

    def solve(self):
        self.prepare_cycles()
        if self.timed_out:
            return

        self.create_aux_graph()
        self.set_lp()
        self.m.set_callbacks(self.callback, lazy=True, cut=True)

        obj_val = -1
        optimal = False
        cuts_correct = False
        while not cuts_correct:
            remaining_time = get_remaining_time(self.start_time)
            if remaining_time > 0.1:
                obj_val, optimal = self.m.solve(remaining_time)
            if not optimal:
                return None
            # cuts_correct = True
            cuts_correct = self.check_cuts()
        return self.get_output(obj_val, optimal), self.__class__

    def set_lp(self):
        n_cycles = len(self.cycles)
        ndds = self.ndds
        self.patients = setdiff(self.vertex_indices, ndds)

        # print("Adding variables")
        self.x = self.m.add_vars(self.edge_tuples)
        self.in_flow = self.m.add_vars(self.vertex_indices)
        self.out_flow = self.m.add_vars(self.vertex_indices)
        self.z = self.m.add_vars(list(range(n_cycles)))

        # print("Adding patient constraints")
        for k in self.patients:
            expr = self.in_flow[k] + -1 * self.out_flow[k]
            if k in self.forbidden_nodes:
                self.m.add_constr_eq(expr, 0)
            else:
                self.m.add_constr_ge(expr, 0)
            cycles_k = self.cycles_part_of[k]
            var_list = [self.z[c] for c in cycles_k]
            var_list.append(self.in_flow[k])
            self.m.add_constr_le(self.m.quick_sum(var_list), 1)

        # print("Adding NDD variables and constraints")
        for n in ndds:
            self.m.add_constr_eq(self.in_flow[n], 0)

        if self.max_chain_length < self.n-1:  # bounded chains
            x_ndds = [None]*len(ndds)
            in_ndds = [None]*len(ndds)
            out_ndds = [None]*len(ndds)
            for n in ndds:
                x_ndds[n] = self.m.add_vars(self.edge_tuples)
                in_ndds[n] = self.m.add_vars(self.vertex_indices)
                out_ndds[n] = self.m.add_vars(self.vertex_indices)

                self.m.add_constr_le(self.m.quick_sum(x_ndds[n].values()), self.max_chain_length)
                for k in self.vertex_indices:
                    self.m.add_constr_eq(self.m.quick_sum([x_ndds[n][(i, k)] for i in self.in_neighbors[k]]) + -1*in_ndds[n][k], 0)
                    self.m.add_constr_eq(self.m.quick_sum([x_ndds[n][(k, j)] for j in self.out_neighbors[k]]) + -1*out_ndds[n][k], 0)
                for p in self.patients:
                    self.m.add_constr_le(out_ndds[n][p] + -1*in_ndds[n][p], 0)
            for e in self.edge_tuples:
                self.m.add_constr_eq(self.m.quick_sum([x_ndds[n][e] for n in ndds]) + -1*self.x[e], 0)
        elif len(ndds) == 0:
            for e in self.edge_tuples:
                self.m.add_constr_eq(self.x[e], 0)

        # print("Setting objective and flow constraints")
        obj_list = [None]*n_cycles
        for c in range(n_cycles):
            obj_list[c] = self.z[c] * self.cycle_weights[c]

        for k in self.vertex_indices:
            expr_in = self.m.quick_sum([self.x[(i, k)] for i in self.in_neighbors[k]]) + -1*self.in_flow[k]
            expr_out = self.m.quick_sum([self.x[(k, j)] for j in self.out_neighbors[k]]) + -1*self.out_flow[k]
            self.m.add_constr_eq(expr_in, 0)
            self.m.add_constr_eq(expr_out, 0)
            for j in self.out_neighbors[k]:
                obj_list.append(self.x[(k, j)] * self.weights[k, j])

        self.m.set_objective_list(obj_list)
        return

    def set_neighbors(self):
        n = self.n

        self.in_neighbors = [None]*n
        self.out_neighbors = [None]*n
        for k in self.vertex_indices:
            self.in_neighbors[k] = self.graph.neighbors(k, mode=ig.IN)
            self.out_neighbors[k] = self.graph.neighbors(k, mode=ig.OUT)

        self.in_neighbors_tuples = [None] * n
        self.out_neighbors_tuples = [None] * n
        for v in range(n):
            self.in_neighbors_tuples[v] = []
            self.out_neighbors_tuples[v] = []
        for e in self.edge_tuples:
            orig = e[0]
            dest = e[1]
            self.in_neighbors_tuples[dest] += [e]
            self.out_neighbors_tuples[orig] += [e]

    def prepare_cycles(self):
        print("Preparing cycles")
        remaining_time = get_remaining_time(self.start_time)
        max_cycle_length = min(self.max_cycle_length, self.n - len(self.ndds))
        cycles = simple_cycles_limited_length(self.graph, max_cycle_length, remaining_time)
        # cycles = simple_cycles_parallel(self.graph, max_cycle_length)
        if cycles is None:
            self.timed_out = True
            return
        num_cycles = len(cycles)
        print("Found %s cycles of length at most %s" % (num_cycles, max_cycle_length))
        w = [0]*num_cycles
        weights = self.weights
        cycles_part_of = [[] for _ in range(self.n)]
        for i, cycle in enumerate(cycles):
            n = len(cycle)
            w[i] = sum(weights[cycle[j], cycle[j+1]] for j in range(n-1))
            w[i] += weights[cycle[n-1], cycle[0]]

            for j in cycle:
                cycles_part_of[j].append(i)
        self.cycles = cycles
        self.cycle_weights = w
        self.cycles_part_of = cycles_part_of

    def create_aux_graph(self):
        cb_graph = self.graph.copy()
        caps = dict()
        n = self.n
        cb_graph.add_vertex()
        new_edges = [None]*len(self.ndds)
        for j, i in enumerate(self.ndds):
            new_edges[j] = (n, i)
            caps[n, i] = 1
            caps[i, n] = math.inf
        cb_graph.add_edges(new_edges)
        self.cb_graph = cb_graph
        self.cb_tuples = [e.tuple for e in cb_graph.es]
        self.added_sets = dict()
        self.caps = caps

    def callback(self, data):
        if len(self.ndds) == 0:
            return
        x_vals = self.m.callback.get_val(self.x, data)
        in_vals = self.m.callback.get_val(self.in_flow, data)
        pos_index = []
        pos_values = []
        for i in self.patients:
            if in_vals[i] > eps:
                pos_index.append(i)
                pos_values.append(in_vals[i])
        edges_sets = self.find_sets(x_vals, pos_index, pos_values)
        for i, edges_set in edges_sets:
            expr = self.m.quick_sum([self.x[e] for e in edges_set]) + -1*self.in_flow[i]
            self.m.callback.add_constr_ge(expr, 0, data)

    def check_cuts(self):
        if len(self.ndds) == 0:
            return True
        x_vals = self.m.get_val(self.x)
        in_vals = self.m.get_val(self.in_flow)
        pos_index = []
        pos_values = []
        for i in self.patients:
            if in_vals[i] > eps:
                pos_index.append(i)
                pos_values.append(in_vals[i])
        edges_sets = self.find_sets(x_vals, pos_index, pos_values)
        for i, edges_set in edges_sets:
            expr = self.m.quick_sum([self.x[e] for e in edges_set]) + -1*self.in_flow[i]
            self.m.add_constr_ge(expr, 0)
        return len(edges_sets) == 0

    def find_sets(self, x_vals, nodes, in_vals):
        n, ndds, graph, cb_graph, caps = self.n, self.ndds, self.graph, self.cb_graph, self.caps
        edges_sets = []
        for i in ndds:
            x_vals[(n, i)] = 1
        edge_tuples = self.cb_tuples
        values = [x_vals[e] for e in edge_tuples]
        cb_graph.es["capacity"] = values
        for v, in_value in zip(nodes, in_vals):
            mincut = cb_graph.mincut(n, v, "capacity")
            flow = mincut.value
            part2 = mincut[1]
            if flow >= in_value-eps:
                continue
            edges_set = []
            for j in setdiff(part2, ndds):
                in_neighbors = self.in_neighbors[j]
                vec = setdiff(in_neighbors, part2)
                for k in vec:
                    edges_set.append((k, j))
            edges_set.sort()  # to avoid repeating constraints
            edges_sets.append((v, edges_set))
        return edges_sets

    def get_output(self, value, optimal):
        match_edges = self.m.get_val(self.x)
        match_cycles = self.m.get_val(self.z)
        graph_cycles = self.cycles
        output = Output(match_edges, match_cycles, graph_cycles, value, optimal)
        return output
