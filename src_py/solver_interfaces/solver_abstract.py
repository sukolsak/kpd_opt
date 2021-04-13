from abc import ABC, abstractmethod


class AbstractSolver(ABC):
    @staticmethod
    @abstractmethod
    def get_solver_name():
        pass

    @staticmethod
    @abstractmethod
    def check_integrality():
        pass

    @abstractmethod
    def add_vars(self, collection: {dict, list}):
        pass

    @abstractmethod
    def set_objective_list(self, obj_list: list):
        pass

    @abstractmethod
    def add_constr_eq(self, expr, rhs: float):
        pass

    @abstractmethod
    def add_constr_le(self, expr, rhs: float):
        pass

    @abstractmethod
    def add_constr_ge(self, expr, rhs: float):
        pass

    @abstractmethod
    def add_sos_constr(self, var_list, weights):
        pass

    @staticmethod
    @abstractmethod
    def quick_sum(var_list):
        pass

    @abstractmethod
    def solve(self, timeout: float):
        pass

    @abstractmethod
    def get_objective_value(self):
        pass

    @staticmethod
    @abstractmethod
    def get_val(collection) -> dict:
        pass

    @abstractmethod
    def set_callbacks(self, function, lazy, cut):
        pass

    @abstractmethod
    def get_n_vars(self):
        pass

    @abstractmethod
    def get_n_constrs(self):
        pass

    @abstractmethod
    def update(self):
        pass


class AbstractCallback(ABC):
    @abstractmethod
    def get_val(self, collection, data):
        pass

    @abstractmethod
    def add_constr_eq(self, expr, rhs, data):
        pass

    @abstractmethod
    def add_constr_le(self, expr, rhs, data):
        pass

    @abstractmethod
    def add_constr_ge(self, expr, rhs, data):
        pass
