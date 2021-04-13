# Solver selection:
from solver_interfaces.solver_glpk import SolverGLPK; solver_instance = SolverGLPK
# from solver_interfaces.solver_cbc import SolverCBC; solver_instance = SolverCBC
# from solver_interfaces.solver_gurobi import SolverGurobi; solver_instance = SolverGurobi
# from solver_interfaces.solver_cplex import SolverCPLEX; solver_instance = SolverCPLEX

# timeout for formulations (in seconds):
timeout = 60
# timeout = int(1e6)

# numerical precision:
eps = 1e-6
