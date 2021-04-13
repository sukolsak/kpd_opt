import time
import igraph as ig

from constants import eps
from formulations.formulation_abstract import Formulation
from utils.graph_utils import find_recipients
from utils.transport import Input, Output
from utils.utils import get_remaining_time, setdiff


class Basic(Formulation):
    def __init__(self, input_data: Input):
        self.solver_instance = input_data.solver_instance
        self.start_time = time.perf_counter()

        self.m = self.solver_instance("basic")
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
        self.set_lp()
        # self.m.set_callbacks(self.callback, lazy=True, cut=False)  # disabling this seems faster for CBC
        return

    def set_neighbors(self):
        self.in_neighbors = [None] * self.n
        self.out_neighbors = [None] * self.n
        for k in self.vertex_indices:
            self.in_neighbors[k] = self.graph.neighbors(k, mode=ig.IN)
            self.out_neighbors[k] = self.graph.neighbors(k, mode=ig.OUT)

        self.in_neighbors_tuples = [None] * self.n
        self.out_neighbors_tuples = [None] * self.n
        for v in range(self.n):
            self.in_neighbors_tuples[v] = []
            self.out_neighbors_tuples[v] = []
        for e in self.edge_tuples:
            orig = e[0]
            dest = e[1]
            self.in_neighbors_tuples[dest] += [e]
            self.out_neighbors_tuples[orig] += [e]

    def set_lp(self):
        # print("Adding variables")
        self.x = self.m.add_vars(self.edge_tuples)
        self.y = self.m.add_vars(self.edge_tuples)
        self.z = self.m.add_vars(self.edge_tuples)

        for e in self.edge_tuples:
            self.m.add_constr_eq(self.x[e] + self.y[e] + -1*self.z[e], 0)

        self.in_flow = self.m.add_vars(self.vertex_indices)
        self.out_flow = self.m.add_vars(self.vertex_indices)

        # print("Adding flow constraints")
        for v in self.vertex_indices:
            in_neighbors_vars = [self.x[e] for e in self.in_neighbors_tuples[v]]
            out_neighbors_vars = [self.x[e] for e in self.out_neighbors_tuples[v]]
            self.m.add_constr_eq(self.m.quick_sum(in_neighbors_vars) + self.in_flow[v] * -1, 0)
            self.m.add_constr_eq(self.m.quick_sum(out_neighbors_vars) + self.out_flow[v] * -1, 0)

        # print("Adding patient constraints")
        ndds = self.ndds
        self.patients = setdiff(self.vertex_indices, ndds)
        for v in self.patients:
            expr = self.out_flow[v] + -1*self.in_flow[v]
            if v in self.forbidden_nodes:
                self.m.add_constr_eq(expr, 0)
            else:
                self.m.add_constr_le(expr, 0)

        # by default chains are unbounded
        obj_terms = [0]*len(self.edge_tuples)
        if self.max_chain_length == 0 or len(ndds) == 0:  # no possible chains, add ad-hoc constraints
            self.cycle_fallback()
            self.set_y_zero()
        elif self.max_chain_length < self.n-1:  # bounded chains, add ad-hoc constraints and new chain variables/constraints
            self.cycle_fallback()
            self.set_ndd_constraints()
            for e in self.edge_tuples:
                edge_vars = [self.x_ndds[v][e] for v in ndds]
                self.m.add_constr_eq(self.m.quick_sum(edge_vars) + -1*self.y[e], 0)
        else:
            self.set_y_zero()

        # print("Setting objective")
        obj_list = [self.z[t] * self.weights[t[0], t[1]] for t in self.edge_tuples if self.weights[t[0], t[1]] != 0]
        self.m.set_objective_list(obj_list)

    def set_y_zero(self):
        for e in self.edge_tuples:
            self.m.add_constr_eq(self.y[e], 0)

    def cycle_fallback(self):
        for p in self.patients:
            self.m.add_constr_eq(self.out_flow[p] + -1*self.in_flow[p], 0)
        for n in self.ndds:
            self.m.add_constr_eq(self.out_flow[n], 0)

    def set_ndd_constraints(self):
        ndds = self.ndds

        x_ndds = [None]*len(ndds)
        in_ndds = [None]*len(ndds)
        out_ndds = [None]*len(ndds)

        for n in ndds:
            x_ndds[n] = self.m.add_vars(self.edge_tuples)
            in_ndds[n] = self.m.add_vars(self.vertex_indices)
            out_ndds[n] = self.m.add_vars(self.vertex_indices)

            for m in ndds:  # break symmetry of flows from NDDs
                self.m.add_constr_eq(in_ndds[n][m], 0)
                if m == n:
                    continue

            self.m.add_constr_le(self.m.quick_sum(x_ndds[n].values()), self.max_chain_length)

            for k in self.vertex_indices:
                self.m.add_constr_eq(self.m.quick_sum([x_ndds[n][(i, k)] for i in self.in_neighbors[k]]) + -1*in_ndds[n][k], 0)
                self.m.add_constr_eq(self.m.quick_sum([x_ndds[n][(k, j)] for j in self.out_neighbors[k]]) + -1*out_ndds[n][k], 0)
            for p in self.patients:
                expr = out_ndds[n][p] + -1*in_ndds[n][p]
                if p in self.forbidden_nodes:
                    self.m.add_constr_eq(expr, 0)
                else:
                    self.m.add_constr_le(expr, 0)

        for k in self.vertex_indices:
            in_flow_ndds = self.m.quick_sum([in_ndds[n][k] for n in ndds])
            out_flow_ndds = self.m.quick_sum([out_ndds[n][k] for n in ndds])
            self.m.add_constr_le(self.in_flow[k] + in_flow_ndds, 1)
            self.m.add_constr_le(self.out_flow[k] + out_flow_ndds, 1)

        self.x_ndds = x_ndds

    def get_match_edges(self, callback=False, data=None):
        if callback:
            z_val = self.m.callback.get_val(self.z, data)
            if z_val is None:
                return None
        else:
            z_val = self.m.get_val(self.z)
        return z_val

    @staticmethod
    def is_integral(var):
        for key, value in var.items():
            if eps < value < 1-eps:
                return False
        return True

    def solve(self):
        found = True
        cycles = []
        obj_val = -1
        optimal = False
        while found:
            self.add_cycle_constrs(cycles)

            optimal = False
            remaining_time = get_remaining_time(self.start_time)
            if remaining_time > 0.1:
                # print("Re solving")
                obj_val, optimal = self.m.solve(remaining_time)
            if not optimal:
                return None
            z_val = self.get_match_edges()
            found, cycles = self.check_cycle_lengths(z_val)
            # print(round(self.m.get_objective_value(), 5), cycles)

        # print("final solution:"); self.print_vars()

        return self.get_output(obj_val, optimal), self.__class__

    def print_vars(self, callback=False, data=None):
        var_dict = {'x': self.x, 'y': self.y, 'z': self.z, 'out:': self.out_flow}
        for key, value in var_dict.items():
            self.print_var(key, value, callback, data)

    def print_var(self, name, var, callback=False, data=None):
        var_dict = self.m.callback.get_val(var, data) if callback else self.m.get_val(var)
        print(name, {key: round(value, 5) for key, value in var_dict.items() if value > eps})

    def get_output(self, value, optimal):
        match_edges = self.get_match_edges()
        match_cycles = None
        graph_cycles = None
        output = Output(match_edges, match_cycles, graph_cycles, value, optimal)
        return output

    def callback(self, data):
        z_val = self.get_match_edges(callback=True, data=data)
        if z_val is None:
            return False

        # print(round(self.m.get_objective_value(), 5))
        # print("current solution:"); self.print_vars(callback=True, data=data)

        if self.solver_instance.check_integrality() and not self.is_integral(z_val):
            return False
        found, cycles = self.check_cycle_lengths(z_val)
        if found:
            self.add_cycle_constrs(cycles, callback=True, data=data)
        return found

    def add_cycle_constrs(self, cycles: list, callback=False, data=None):
        for element in cycles:
            cycle_length = element[0]
            cycle = element[1]

            cycle_vars = [self.z[e] for e in cycle]
            expr = self.m.quick_sum(cycle_vars)
            if callback:
                self.m.callback.add_constr_le(expr, cycle_length-1, data)
            else:
                self.m.add_constr_le(expr, cycle_length-1)

    def check_cycle_lengths(self, x_val, early_stop=False):
        n = self.n
        flag = [False]*n
        cycles = []
        recipients = find_recipients(self.graph, x_val)
        for base in range(n):
            if flag[base]:
                continue
            current_vertex = base
            cycle_length = 0
            cycle = []
            while True:
                flag[current_vertex] = True
                next_vertex = recipients[current_vertex]
                if next_vertex == -1:  # no recipient, end of a chain
                    break
                flag[next_vertex] = True
                cycle += [(current_vertex, next_vertex)]
                cycle_length += 1
                if next_vertex == base:  # found the cycle
                    if cycle_length > self.max_cycle_length:
                        cycles += [(cycle_length, cycle)]
                        if early_stop:
                            return True, cycles
                    break
                current_vertex = next_vertex
        return len(cycles) > 0, cycles
