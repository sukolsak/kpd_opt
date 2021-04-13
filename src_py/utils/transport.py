from typing import List, Dict, Tuple
import igraph as ig
import numpy as np


class Input:
    def __init__(
            self,
            graph: ig.Graph,
            weights: np.ndarray,
            ndds: List[int],
            cycle_length: int,
            chain_length: int,
            forbidden_nodes: List[int] = None
            ):
        self.graph = graph
        self.weights = weights
        self.ndds = ndds
        self.cycle_length = cycle_length
        self.chain_length = chain_length
        self.solver_instance = None
        self.start_time = None

        self.forbidden_nodes = forbidden_nodes


class Output:
    def __init__(
            self,
            match_edges: Dict[Tuple[int], float],
            match_cycles: {Dict[int, float], None},
            graph_cycles: {List[Tuple[int]], None},
            value: float,
            optimal: bool,
            ):
        self.match_edges = match_edges
        self.match_cycles = match_cycles
        self.graph_cycles = graph_cycles
        self.value = value
        self.optimal = optimal
        self.time = -1
