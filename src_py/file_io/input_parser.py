import math


def parse_input(dirname:str):
    with open(dirname) as f:
        lines = f.readlines()
    assert lines[1] == "problemData"
    chain_length = lines[2].split(',')
    assert chain_length[1] == "maxChainLength"
    if chain_length[2] == "Infinity":
        max_chain_length = math.inf
    else:
        max_chain_length = int(chain_length[2])
    return max_chain_length
