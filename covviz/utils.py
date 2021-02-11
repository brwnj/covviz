import gzip

import numpy as np


def gzopen(f):
    if f.endswith(".gz"):
        return gzip.open(f, "rt")
    else:
        return open(f)


def compare_array(a, b):
    len_a = len(a)
    len_b = len(b)
    if len_a > len_b:
        b = np.pad(b, (0, len_a - len_b), mode="constant")
    elif len_b > len_a:
        a = np.pad(a, (0, len_b - len_a), mode="constant")
    try:
        return a[: np.nonzero(a != b)[0][0]]
    except IndexError:
        # case: a == b
        return a


def find_common_start(items):
    if len(items) == 1:
        return items

    res = compare_array(np.array(items[0]), np.array(items[1]))
    if len(items) > 2:
        for i in range(2, len(items)):
            res = compare_array(res, np.array(items[i]))
    return res.tolist()


def optimize_coords(data):
    """
    adds list to incoming dict (data) under key "shared_coords"
    per chrom, shortens coords to remove redundant, shared starting elements
    """
    coords = []
    for chrom in data["chromosomes"]:
        coords.append(data[chrom]["coords"])
    # find common starting elements
    data["shared_coords"] = find_common_start(coords)
    len_coords = len(data["shared_coords"])
    # remove common starting elements
    for chrom in data["chromosomes"]:
        data[chrom]["coords"] = data[chrom]["coords"][len_coords:]
    return data
