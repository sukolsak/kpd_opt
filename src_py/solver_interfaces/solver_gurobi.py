import math

from gurobipy import *

from solver_interfaces.solver_abstract import AbstractSolver, AbstractCallback

setParam('OutputFlag', 0)
setParam('LazyConstraints', 1)


class SolverGurobi(AbstractSolver):
    solver_name = 'Gurobi'

    def __init__(self, name):
        self.m = Model(name)
        self.callback = None

    @staticmethod
    def get_solver_name():
        return SolverGurobi.solver_name

    @staticmethod
    def check_integrality():
        return False

    def add_vars(self, collection: {dict, list}):
        keys = None
        if isinstance(collection, dict):
            keys = collection.keys()
        elif isinstance(collection, list):
            keys = collection
        new_vars = self.m.addVars(keys, vtype=GRB.BINARY)
        return new_vars

    def set_objective_list(self, obj_list: list):
        self.m.setObjective(quicksum(obj_list), sense=GRB.MAXIMIZE)

    def add_constr_eq(self, expr, rhs: float):
        self.m.addConstr(expr == rhs)

    def add_constr_le(self, expr, rhs: float):
        self.m.addConstr(expr <= rhs)

    def add_constr_ge(self, expr, rhs: float):
        self.m.addConstr(expr >= rhs)

    def add_sos_constr(self, var_list, weights):
        self.m.addSOS(GRB.SOS_TYPE1, var_list, weights)

    @staticmethod
    def quick_sum(var_list):
        return quicksum(var_list)

    def solve(self, timeout=math.inf):
        def callback_function(model, where):
            data = [model, where]
            if self.callback is None or where not in self.callback.reasons_dict.keys():
                return
            violated = self.callback.function(data)
            return
        self.m.setParam('TimeLimit', timeout)
        self.m.optimize(callback_function)
        optimal = self.m.Status == GRB.OPTIMAL
        return self.m.objVal, optimal

    def get_objective_value(self):
        return self.m.objVal

    @staticmethod
    def get_val(collection: tupledict) -> dict:
        vars_val = {key: v.X for key, v in collection.items()}
        return vars_val

    def set_callbacks(self, function, lazy=True, cut=True):
        reasons_dict = dict()
        if lazy:
            reasons_dict[GRB.Callback.MIPSOL] = "GRB.Callback.MIPSOL"
        if cut:
            reasons_dict[GRB.Callback.MIPNODE] = "GRB.Callback.MIPNODE"
        self.callback = CallbackGurobi(self, function, reasons_dict)
        return

    def get_n_vars(self):
        return self.m.NumVars

    def get_n_constrs(self):
        return self.m.NumConstrs

    def update(self):
        self.m.update()


class CallbackGurobi(AbstractCallback):
    def __init__(self, model, function, reasons_dict):
        self.model = model
        self.function = function
        self.reasons_dict = reasons_dict
        return

    def get_val(self, collection: tupledict, data) -> dict:
        model, where = data
        indices = list(collection.keys())
        vars_list = [collection[key] for key in indices]
        # TODO: test extensively if this is the correct flow for callbacks
        if where == GRB.Callback.MIPNODE:
            status = model.cbGet(GRB.Callback.MIPNODE_STATUS)
            values_list = model.cbGetNodeRel(vars_list)
            # if status == GRB.OPTIMAL:
            #     values_list = model.cbGetNodeRel(vars_list)
            # else:
            #     values_list = model.cbGetNodeRel(vars_list)
            #     return None
        elif where == GRB.Callback.MIPSOL:
            values_list = model.cbGetSolution(vars_list)
        else:
            values_list = model.cbGetSolution(vars_list)
        return {key: values_list[i] for i, key in enumerate(indices)}

    @staticmethod
    def add_constr(expr, data):
        model, where = data
        if where == GRB.Callback.MIPNODE:
            model.cbCut(expr)
        elif where == GRB.Callback.MIPSOL:
            model.cbLazy(expr)

    def add_constr_le(self, expr, rhs, data):
        self.add_constr(expr <= rhs, data)

    def add_constr_eq(self, expr, rhs: float, data):
        self.add_constr(expr == rhs, data)

    def add_constr_ge(self, expr, rhs: float, data):
        self.add_constr(expr >= rhs, data)
