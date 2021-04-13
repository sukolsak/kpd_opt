import multiprocessing
from typing import Union, Tuple

from cplex.callbacks import UserCutCallback, LazyConstraintCallback
from docplex.mp.callbacks.cb_mixin import ConstraintCallbackMixin
from docplex.mp.linear import LinearExpr
from docplex.mp.model import Model
from docplex.mp.solution import SolveSolution

from solver_interfaces.solver_abstract import AbstractSolver, AbstractCallback


class SolverCPLEX(AbstractSolver):
    solver_name = 'CPLEX'

    def __init__(self, name):
        model = Model(name)
        model.parameters.threads = multiprocessing.cpu_count()

        self.m = model
        self.callback = None

    @staticmethod
    def get_solver_name():
        return SolverCPLEX.solver_name

    @staticmethod
    def check_integrality():
        return False

    def add_vars(self, collection: {dict, list}):
        keys = None
        if isinstance(collection, dict):
            keys = collection.keys()
        elif isinstance(collection, list):
            keys = collection
        new_vars = self.m.binary_var_dict(keys)
        return new_vars

    def set_objective_list(self, obj_list: list):
        self.m.set_objective("max", self.quick_sum(obj_list))

    def add_constr_eq(self, expr, rhs: float):
        self.m.add_constraint(expr == rhs)

    def add_constr_le(self, expr, rhs: float):
        self.m.add_constraint(expr <= rhs)

    def add_constr_ge(self, expr, rhs: float):
        self.m.add_constraint(expr >= rhs)

    def add_sos_constr(self, var_list, weights):
        vars_sorted = [var for _, var in sorted(zip(weights, var_list))]
        self.m.add_sos1(vars_sorted)

    def quick_sum(self, var_list):
        return self.m.sum(var_list)

    def solve(self, timeout: float):
        if timeout < 1:
            timeout = 1
        self.m.set_time_limit(int(timeout))
        self.m.solve()
        optimal = self.m.get_solve_status() is not None
        return self.m.objective_value, optimal

    def get_objective_value(self):
        return self.m.objective_value

    @staticmethod
    def get_val(collection) -> dict:
        vars_val = {key: v.solution_value for key, v in collection.items()}
        return vars_val

    def set_callbacks(self, function, lazy, cut):
        self.callback = CallbackCPLEX(self, function)
        test_cb = None
        if lazy:
            lazy_cb = self.m.register_callback(CallbackCPLEXLazy)
            lazy_cb.general_cb = self.callback
        if cut:
            cut_cb = self.m.register_callback(CallbackCPLEXCut)
            cut_cb.general_cb = self.callback
        return

    def get_n_vars(self):
        return self.m.number_of_variables

    def get_n_constrs(self):
        return self.m.number_of_constraints

    def update(self):
        return


class CallbackCPLEXCut(ConstraintCallbackMixin, UserCutCallback):
    def __init__(self, env):
        self.general_cb = None
        UserCutCallback.__init__(self, env)
        ConstraintCallbackMixin.__init__(self)

    def __call__(self):
        self.general_cb.execute(self)


class CallbackCPLEXLazy(ConstraintCallbackMixin, LazyConstraintCallback):
    def __init__(self, env):
        self.general_cb = None
        LazyConstraintCallback.__init__(self, env)
        ConstraintCallbackMixin.__init__(self)

    def __call__(self):
        self.general_cb.execute(self)


CPLEXCallback = Union[CallbackCPLEXCut, CallbackCPLEXLazy]
Data = Tuple[CPLEXCallback, SolveSolution]


class CallbackCPLEX(AbstractCallback):
    def __init__(self, model: SolverCPLEX, function):
        self.model = model
        self.function = function

    def execute(self, cplex_cb: CPLEXCallback):
        sol = cplex_cb.make_solution()
        data = [cplex_cb, sol]
        self.function(data)

    def get_val(self, collection, data: Data):
        sol = data[1]
        return sol.get_value_dict(collection)

    @staticmethod
    def add_constr(expr: LinearExpr, sense, rhs, data: Data):
        cplex_cb = data[0]
        cplex_cb.add(expr, sense, rhs)
        '''
        if sense == "L":
            self.model.m.add_constraint(expr <= rhs)
        elif sense == "G":
            self.model.m.add_constraint(expr >= rhs)
        elif sense == "E":
            self.model.m.add_constraint(expr == rhs)
        # '''

    def add_constr_le(self, expr, rhs: float, data):
        self.add_constr(expr, "L", rhs, data)

    def add_constr_eq(self, expr, rhs: float, data):
        self.add_constr(expr, "E", rhs, data)

    def add_constr_ge(self, expr, rhs: float, data):
        self.add_constr(expr, "G", rhs, data)
