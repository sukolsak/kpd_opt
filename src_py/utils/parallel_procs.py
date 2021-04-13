from multiprocessing import Pool


def from_list(fun, input_list: list, chunksize: int = None, processes: int = None) -> list:
    input_list = [[fun] + e for e in input_list]
    pool = Pool(processes=processes)
    try:
        result_list = pool.map(function_master, input_list, chunksize=chunksize)
    finally:
        pool.close()
        pool.join()
    return result_list


def from_list_series(fun, input_list: list) -> list:
    input_list = [[fun, False] + e for e in input_list]
    return list(map(function_master, input_list))


def function_master(work_data):
    fun = work_data[0]
    input_tuple = tuple(work_data[1:])
    result = fun(*input_tuple)
    return result
