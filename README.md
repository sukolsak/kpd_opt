# kpd_opt

**kpd_opt** is a tool for finding optimal matching between donors and recipients in kidney paired donation. It is an implementation of the algorithms described in [the paper "Finding long chains in kidney exchange using the traveling salesman problem"](https://web.stanford.edu/~iashlagi/papers/pnasChain.pdf) by Anderson et al. (2015). It takes as input a compatibility graph where each node represents an incompatible donor-recipient pair or a non-directed donor and each edge represents donor-recipient compability. It outputs cycles and chains of exchanges.

## Getting started, Julia version
1. Install Julia and GLPK
2. Install the Julia packages to deal with graphs:
```
pkg> add LightGraphs
pkg> add LightGraphsFlows
pkg> add LightGraphsMatching
```
3. Install the Julia packages to interface with GLPK:
```
pkg> add GLPK
pkg> add GLPKMathProgInterface
```
4. Install the latest version of JuMP (from github) to solve with callbacks:
```
pkg> add https://github.com/JuliaOpt/JuMP.jl
```
5. Run from terminal:
```
julia main.jl <input path> <output path> <formulation>
```
The formulation may be

- `basic` for the recursive formulation
- `pctsp` to use the PC-TSP formulation
- `parallel` to run both in parallel processes

See the example input and output files in [/examples](examples).

## Getting started, Python version
1. Install Python 3 and GLPK
2. Install the following packages:
 - [numpy](https://github.com/numpy/numpy): Scientific library, used to work with the weight matrices.
 - [lap](https://github.com/gatagat/lap): Fast C++ implementation of the Jonker-Volgenant algorithm to solve the assignment problem, used by the fallback formulation.
 - [python-igraph](https://igraph.org/python/): Fast C/C++ library to work with graphs, used to work with in/out-neighbors, solve the maximum matching problem and the minimum cut subproblems.
 - [glpk](http://tfinley.net/software/pyglpk/): Python interface to GLPK solver, used to solve the MIPs with callbacks.
```
pip install numpy
pip install lap
pip install python-igraph
pip install glpk
```
3. Run from terminal:
```
python main.py <input path> <output path> <formulation>
```
The formulation may be

- `basic` for the recursive formulation
- `pctsp` to use the PC-TSP formulation
- `fallback` to use the combinatorial matching formulation if applicable (unbounded cycles and chains, or 2-cycles and no chains)
- `parallel` to run all in parallel processes

You can select the solver in [constants.py](src_py/constants.py). GLPK, Cbc, Gurobi, and CPLEX are supported.

See the example input and output files in [/examples](examples).

## Contributors

Itai Ashlagi, Lilia Chang, Sukolsak Sakshuwong, Felipe Subiabre
