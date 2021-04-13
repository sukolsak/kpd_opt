import multiprocessing as mp
import time

from constants import timeout
from formulations.basic import Basic
from formulations.fallback import Fallback
from formulations.formulation_abstract import Formulation
from formulations.intermediate import Intermediate

formulations_list = [Fallback, Basic, Intermediate]
n = len(formulations_list)

dt = 0.01


class Parallel(Formulation):
    def __init__(self, *args):
        self.args = args
        self.procs = [None]*n
        self.start_time = time.perf_counter()

    def solve(self):
        sol_queue = mp.Queue()
        finished_queue = mp.Queue(maxsize=n)
        procs = self.procs
        for i, formulation in enumerate(formulations_list):
            procs[i] = mp.Process(target=spawn, args=(formulation, self.args, sol_queue, finished_queue))
            procs[i].start()

        result = None
        tf = self.start_time
        while tf - self.start_time < timeout:
            if not sol_queue.empty():
                result = sol_queue.get()
                break
            if finished_queue.full():
                break
            time.sleep(dt)
            tf = time.perf_counter()
        self.terminate_procs()
        return result

    def terminate_procs(self):
        for p in self.procs:
            p.terminate()


def spawn(formulation, args, sol_queue: mp.Queue, finished_queue: mp.Queue):
    name = formulation.__name__
    # print(name, "started")
    problem = formulation(*args)
    # print(name, "optimizing")
    result = problem.solve()
    if result is not None:
        sol_queue.put(result)
        print(name, "finished")
    finished_queue.put(formulation)
