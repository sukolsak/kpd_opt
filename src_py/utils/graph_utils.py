import math
import time
import igraph as ig

from utils.parallel_procs import from_list
# from parallel_procs import from_list


def find_cycles(g: ig.Graph, k, time_limit=math.inf):  # own implementation
    t0 = time.perf_counter()
    path_dict = dict()
    cycles = []
    in_neighs = [g.neighbors(v, mode=ig.IN) for v in g.vs]
    out_neighs = [g.neighbors(v, mode=ig.OUT) for v in g.vs]
    for length in range(k):
        for u in g.vs:
            tf = time.perf_counter()
            if tf-t0 > time_limit:
                return None
            u_index = u.index
            paths = find_simple_paths_node(g, u_index, length, path_dict, out_neighs)
            in_neighbors = in_neighs[u_index]
            cycles_u = [list(path) for path in paths if path[-1] in in_neighbors]
            for cycle in cycles_u:
                add = True
                for v in cycle[1:]:
                    if v >= u_index:
                        add = False
                        break
                if add:
                    cycles += [cycle]
    return sorted(cycles)


def find_simple_paths_node(g: ig.Graph, u, k, path_dict=None, out_neigh=None):
    paths_uk = path_dict.get((u, k), None)
    if paths_uk is not None:
        return paths_uk
    if k == 0:
        paths_uk = [[u]]
    else:
        paths_uk = []
        for v in out_neigh[u]:
            paths_vk = path_dict.get((v, k-1), None)
            if paths_vk is None:
                paths_vk = []
                paths_vk += find_simple_paths_node(g, v, k-1, path_dict, out_neigh)
            paths_uk += [[u]+p for p in paths_vk if u not in p]
    path_dict[(u, k)] = paths_uk
    return paths_uk


def simple_cycles_limited_length(graph: ig.Graph, k, remaining_time=math.inf):  # adapted julia implementation
    max_time = time.perf_counter()+remaining_time
    cycles = []
    if k < 1:
        return cycles
    cycle = [None]*k
    out_neighs = [graph.neighbors(v, ig.OUT) for v in graph.vs]
    for v in graph.vs:
        if time.perf_counter() > max_time:
            return None
        cycle[0] = v.index
        status = simple_cycles_limited_length_node(graph, k, cycles, cycle, 0, out_neighs)
        if status == -1:
            return None
    return sorted(cycles)


def simple_cycles_limited_length_node(graph: ig.Graph, k, cycles, cycle, i, out_neighs):
    neighbors = out_neighs[cycle[i]]
    repeated = False
    for v in neighbors:
        if v == cycle[0] and not repeated:
            cycles.append(cycle[:i+1])
            repeated = True
        elif i < k-1 and v > cycle[0] and not repeated_vertex(v, cycle, 1, i):
            cycle[i+1] = v
            simple_cycles_limited_length_node(graph, k, cycles, cycle, i+1, out_neighs)
    return 0


def simple_cycles_parallel(graph: ig.Graph, k):
    if k < 1:
        return []
    out_neighs = [graph.neighbors(v, ig.OUT) for v in graph.vs]
    input_list = [None]*graph.vcount()
    for v in graph.vs:
        cycle = [None]*k
        cycle[0] = v.index
        input_list[v.index] = [graph, k, [], cycle, 0, out_neighs]
    cycles_node = from_list(simple_cycles_parallel_node, input_list)
    cycles = list()
    for e in cycles_node:
        cycles.extend(e)
    return sorted(cycles)


def simple_cycles_parallel_node(graph: ig.Graph, k, cycles, cycle, i, out_neighs):
    neighbors = out_neighs[cycle[i]]
    repeated = False
    for v in neighbors:
        if v == cycle[0] and not repeated:
            cycles.append(cycle[:i+1])
            repeated = True
        elif i < k-1 and v > cycle[0] and not repeated_vertex(v, cycle, 1, i):
            cycle[i+1] = v
            simple_cycles_parallel_node(graph, k, cycles, cycle, i+1, out_neighs)
    return cycles


def repeated_vertex(v, cycle, n1, n2):
    for k in range(n1, n2+1):
        if cycle[k] == v:
            return True
    return False


def find_recipients(graph: ig.Graph, x_val: dict):
    recipients = [-1]*graph.vcount()
    for e, value in x_val.items():
        i = e[0]
        j = e[1]
        if value > 1/2:
            recipients[i] = j
    return recipients


def main():
    n = 154
    m = int(n*n*0.14)
    # m = int(5e3)
    k = 4
    g = ig.Graph.Erdos_Renyi(n, m=m, directed=True, loops=False)

    implementations = [simple_cycles_parallel, simple_cycles_limited_length, find_cycles]

    for f in implementations:
        print(f.__name__)
        t0 = time.perf_counter()
        cycles = f(g, k)
        tf = time.perf_counter()
        print(round(tf-t0, 3), 'sec')
        print(len(cycles), 'cycles')
        print()


if __name__ == '__main__':
    main()
