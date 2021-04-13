import time

from constants import timeout


def get_remaining_time(start_time):
    elapsed_time = time.perf_counter() - start_time
    return max(timeout - elapsed_time, 0)


def setdiff(a, b):
    return list(set(a).difference(b))
