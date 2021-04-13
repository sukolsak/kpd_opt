import math
import sys
import time
from copy import copy

import glpk
from solver_interfaces.solver_abstract import AbstractSolver, AbstractCallback

ignore_msgs = ['Long-step dual simplex will be used']


def terminal_hook(*args):
    for arg in args:
        arg = arg.strip()
        if arg in ignore_msgs:
            continue
        print("GLPK:", arg)


glpk.env.term_hook = terminal_hook


class Variable:
    index = 0

    def __init__(self, var: glpk.Bar):
        self.var = var
        self.index_model = -1
        self.model = None

        self.index = Variable.index
        Variable.index += 1

    def __str__(self):
        name = self.var.name
        if name is None:
            name = "var" + str(self.index)
        return name

    def to_expression(self):
        expr = Expression()
        expr.add_var(self, 1)
        return expr

    def __mul__(self, other: float):
        expr = Expression()
        expr.add_var(self, other)
        return expr

    def __add__(self, other):
        if isinstance(other, int):
            assert other == 0
            return self
        if isinstance(other, Variable):
            new = Expression()
            new.add_var(self, 1)
            new.add_var(other, 1)
        elif isinstance(other, Expression):
            new = copy(other)
            new.add_var(self, 1)
        else:
            sys.exit("Error: can't add variable to object of type %s" % type(other))
        return new

    def __radd__(self, other):
        return self.__add__(other)

    def __rmul__(self, other):
        return self.__mul__(other)


class Expression:
    def __init__(self):
        self.vars = []
        self.coeffs = []
        self.indices = []

    def __len__(self):
        return len(self.vars)

    def __iter__(self):
        return zip(self.vars, self.coeffs)

    def __str__(self):
        return repr(self)

    def __copy__(self):
        new = Expression()
        n = len(self)
        new.coeffs = [None]*n
        new.vars = [None]*n
        new.indices = [None]*n
        for i in range(n):
            new.coeffs[i] = self.coeffs[i]
            new.vars[i] = self.vars[i]
            new.indices[i] = self.indices[i]
        return new

    def __repr__(self):
        text = ""
        i = 0
        for v, coeff in self:
            if i > 0:
                text += " + "
            if coeff == 1:
                text += str(v)
            else:
                text += str(coeff) + "*" + str(v)
            i += 1
        return text

    def add_var(self, var: Variable, coeff: float):
        try:
            i = self.indices.index(var.index)
            self.coeffs[i] += coeff
        except ValueError:
            self.vars += [var]
            self.coeffs += [coeff]
            self.indices += [var.index]

    def __add__(self, other):
        new = copy(self)
        if isinstance(other, int):
            assert other == 0
            return new
        if isinstance(other, Expression):
            i = 0
            for var, coeff in other:
                new.add_var(var, coeff)
                i += 1
        elif isinstance(other, Variable):
            new.add_var(other, 1)
        else:
            sys.exit("Error: can't add expression to object of type %s" % type(other))
        return new

    def __mul__(self, other: float):
        new = copy(self)
        for i in range(len(self)):
            new.coeffs[i] *= other
        return new

    def __radd__(self, other):
        return self.__add__(other)

    def __rmul__(self, other):
        return self.__mul__(other)


class SolverGLPK(AbstractSolver):
    solver_name = 'GLPK'

    @staticmethod
    def get_solver_name():
        return SolverGLPK.solver_name

    @staticmethod
    def check_integrality():
        return True

    def __init__(self, name: str):
        m = glpk.LPX()
        self.m = m
        self.m.name = name
        self.vars = []
        self.callback = None
        self.timeout = math.inf

    def add_vars(self, collection: {dict, list}):
        new_vars = dict()
        if isinstance(collection, dict):
            for key in collection.keys():
                new_vars[key] = self.add_var()
            return new_vars
        elif isinstance(collection, list):
            for v in collection:
                new_vars[v] = self.add_var()
            return new_vars
        return None

    def add_var(self, name: str = None):
        m = self.m
        m.cols.add(1)
        var = m.cols[-1]
        var.kind = bool
        if name is not None:
            var.name = name
        var_ext = Variable(var)
        var_ext.index_model = len(self)-1
        var_ext.model = m
        self.vars += [var_ext]
        return var_ext

    def __len__(self):
        return len(self.m.cols)

    def check_expr(self, expr: Expression):
        m = self.m
        assert all([var.model == m for var, coeff in expr])

    def set_objective_list(self, obj_list: list):
        n = len(self)
        row = [0]*n
        for e in obj_list:
            assert isinstance(e, Expression)
            self.check_expr(e)
            var = e.vars[0]
            coeff = e.coeffs[0]
            row[var.index_model] += coeff

        m = self.m
        m.obj[:] = row
        m.obj.maximize = True

    def set_objective(self, expr: {Variable, Expression}):
        if isinstance(expr, Variable):
            expr = expr.to_expression()
        self.check_expr(expr)
        m = self.m
        n = len(self)
        row = [0]*n
        for var, coeff in expr:
            row[var.index_model] = coeff
        m.obj[:] = row
        m.obj.maximize = True

    def add_constr(self, expr: {Variable, Expression}, rhs_low, rhs_high):
        if isinstance(expr, Variable):
            expr = expr.to_expression()
        self.check_expr(expr)
        m = self.m
        m.rows.add(1)
        row = m.rows[-1]
        row.bounds = rhs_low, rhs_high
        row.matrix = [(var.index_model, coeff) for var, coeff in expr]

    def add_constr_eq(self, expr: {Variable, Expression}, rhs: float):
        self.add_constr(expr, rhs, rhs)

    def add_constr_le(self, expr: {Variable, Expression}, rhs: float):
        self.add_constr(expr, None, rhs)

    def add_constr_ge(self, expr: {Variable, Expression}, rhs: float):
        self.add_constr(expr, rhs, None)

    def add_sos_constr(self, var_list, weights):  # GLPK has no special treatment for SOS constraints
        self.add_constr_le(var_list, 1)

    def set_callbacks(self, function, lazy=True, cut=True):
        reasons = []
        if lazy:
            reasons += ["rowgen"]
        if cut:
            reasons += ["cutgen"]
        self.callback = CallbackGLPK(self, function, reasons)

    def solve(self, timeout=math.inf):
        t0 = time.perf_counter()
        value = -1
        msg_lev = glpk.LPX.MSG_OFF
        # msg_lev = glpk.LPX.MSG_ALL
        self.m.simplex(
            msg_lev=msg_lev, tm_lim=int(timeout*1000),
            presolve=True
        )
        if self.m.status != 'opt':
            return value, False
        t1 = time.perf_counter()
        remaining_time = max(timeout-(t1-t0), 1)
        self.m.integer(
            callback=self.callback, msg_lev=msg_lev, tm_lim=int(remaining_time*1000),
            gmi_cuts=True,
            mir_cuts=True,
            # pp_tech=glpk.LPX.PP_ROOT,
        )
        optimal = self.m.status == 'opt'
        if optimal:
            value = self.m.obj.value
            for v in self.vars:
                v.value = v.var.primal
        return value, optimal

    def get_objective_value(self):
        return self.m.obj.value

    def set_timeout(self, timeout):
        self.timeout = timeout

    @staticmethod
    def quick_sum(var_list: list):
        expr = Expression()
        n = len(var_list)
        expr.vars = var_list
        expr.coeffs = [1]*n
        expr.indices = [var.index for var in var_list]
        return expr

    @staticmethod
    def get_val(collection) -> dict:
        vars_val = {key: v.var.primal for key, v in collection.items()}
        return vars_val

    def get_n_vars(self):
        return len(self.m.cols)

    def get_n_constrs(self):
        return len(self.m.rows)

    def update(self):
        return


class CallbackGLPK(AbstractCallback):
    def __init__(self, model: SolverGLPK, function, reasons):
        self.model = model
        self.function = function
        self.reasons = reasons
        return

    def general(self, tree: glpk.Tree):
        self.function(tree)
        return

    def default(self, tree: glpk.Tree):
        if tree.reason in self.reasons:
            self.general(tree)
        return

    def add_constr(self, expr: {Variable, Expression}, rhs_low, rhs_high, tree):
        self.model.add_constr(expr, rhs_low, rhs_high)
        m = tree.lp
        m.rows.add(1)
        row = m.rows[-1]
        row.bounds = rhs_low, rhs_high
        row.matrix = [(var.index_model, coeff) for var, coeff in expr]

    def add_constr_eq(self, expr: {Variable, Expression}, rhs: float, data):
        self.add_constr(expr, rhs, rhs, data)

    def add_constr_le(self, expr: {Variable, Expression}, rhs: float, data):
        self.add_constr(expr, None, rhs, data)

    def add_constr_ge(self, expr: {Variable, Expression}, rhs: float, data):
        self.add_constr(expr, rhs, None, data)

    def get_val(self, collection, data) -> dict:
        return self.model.get_val(collection)
