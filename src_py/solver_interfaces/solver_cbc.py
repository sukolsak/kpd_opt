import os
import sys

from mip import *
from solver_interfaces.solver_abstract import AbstractSolver, AbstractCallback


class SolverCBC(AbstractSolver):
    solver_name = 'CBC'

    def __init__(self, name):
        model = Model(name, sense=MAXIMIZE, solver_name=CBC)
        # model = Model(name, sense=MAXIMIZE, solver_name=GUROBI)
        model.verbose = 0

        # model.preprocess = 1
        # model.clique = 0
        # model.emphasis = 0
        # model.cuts = -1

        self.m = model
        self.callback = None
        self.errored = False

    @staticmethod
    def get_solver_name():
        return SolverCBC.solver_name

    @staticmethod
    def check_integrality():
        return False

    def add_vars(self, collection: {dict, list}):
        keys = None
        if isinstance(collection, dict):
            keys = collection.keys()
        elif isinstance(collection, list):
            keys = collection
        new_vars = {key: self.m.add_var(var_type=BINARY) for key in keys}
        return new_vars

    def set_objective_list(self, obj_list: list):
        self.obj_list = obj_list
        self.m.objective = maximize(xsum(obj_list))

    def add_constr_eq(self, expr, rhs: float):
        self.m.add_constr(expr == rhs)

    def add_constr_le(self, expr, rhs: float):
        self.m.add_constr(expr <= rhs)

    def add_constr_ge(self, expr, rhs: float):
        self.m.add_constr(expr >= rhs)

    def add_sos_constr(self, var_list, weights):
        sos = list(zip(var_list, weights))
        self.m.add_sos(sos, 1)

    @staticmethod
    def quick_sum(var_list):
        return xsum(var_list)

    def solve(self, timeout: float):
        status = self.m.optimize(max_seconds=timeout)
        optimal = (status == OptimizationStatus.OPTIMAL) and not self.errored
        return self.m.objective_value, optimal

    def get_objective_value(self):
        return self.m.objective_value

    @staticmethod
    def get_val(collection) -> dict:
        vars_val = {key: v.x for key, v in collection.items()}
        return vars_val

    def set_callbacks(self, function, lazy, cut):
        self.callback = CallbackCBC(self, function)
        if lazy:
            self.m.lazy_constrs_generator = self.callback
        if cut:
            self.m.cuts_generator = self.callback
        return

    def get_n_vars(self):
        return self.m.num_cols

    def get_n_constrs(self):
        return self.m.num_rows

    def update(self):
        return


class CallbackCBC(AbstractCallback, ConstrsGenerator):
    def __init__(self, model: SolverCBC, function):
        self.model = model
        self.function = function
        self.cut_pool = CutPool()
        return

    def generate_constrs(self, data: Model, depth: int = 0, npass: int = 0):
        # print(self.get_obj_val(data))
        self.function(data)

    def add_constr(self, expr: LinExpr, data: Model):
        # added = self.cut_pool.add(expr)
        # if not added:
        #     return
        sense = expr.sense
        n = len(expr.expr)
        variables = [None]*n
        coeffs = [0]*n
        for i, (key, value) in enumerate(expr.expr.items()):
            variables[i] = key
            coeffs[i] = value
        new_vars = data.translate(variables)
        new_expr = LinExpr(new_vars, coeffs, expr.const, sense)
        data.add_cut(new_expr)

    def add_constr_le(self, expr, rhs: float, data):
        self.add_constr(expr <= rhs, data)

    def add_constr_eq(self, expr, rhs: float, data):
        self.add_constr(expr == rhs, data)

    def add_constr_ge(self, expr, rhs: float, data):
        self.add_constr(expr >= rhs, data)

    def get_val(self, collection: dict, data: Model):
        # TODO: fix when Python-MIP gets fixed
        collection_new = dict()
        error = False
        for key, value in collection.items():
            new_value = data.translate(value)
            if new_value is None:
                error = True
                new_value = -1
            collection_new[key] = new_value
        if error:
            print("CBC internal error: variable not found in preprocessed model")
            error_var = data.vars[0]
            data.add_constr(error_var >= 1)
            data.add_constr(error_var <= 0)
            self.model.errored = True
            return collection_new
        return self.model.get_val(collection_new)

    def get_obj_val(self, data):
        obj_list = self.model.obj_list
        obj_coeffs = [None]*len(obj_list)
        obj_vars = dict()
        for i, expr in enumerate(obj_list):
            for key, value in expr.expr.items():
                obj_vars[i] = key
                obj_coeffs[i] = value
        values = self.get_val(obj_vars, data)
        obj_val = 0
        for i, value in values.items():
            obj_val += obj_coeffs[i]*value
        return obj_val
