import gzip

# 3 sets of 10 colors
COLORS = [
    "#1f77b4",  # plotly and d3
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#7CB5EC",  # highcharts
    "#434348",
    "#90ED7D",
    "#F7A35C",
    "#8085E9",
    "#F15C80",
    "#E4D354",
    "#2B908F",
    "#F45B5B",
    "#91E8E1",
    "#4E79A7",  # tableau
    "#F28E2C",
    "#E15759",
    "#76B7B2",
    "#59A14F",
    "#EDC949",
    "#AF7AA1",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AB",
]


def gzopen(f):
    if f.endswith(".gz"):
        return gzip.open(f, "rt")
    else:
        return open(f)


def merge_intervals(intervals):
    sorted_intervals = sorted(intervals, key=lambda i: i[0])
    merged = list()
    for higher in sorted_intervals:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # merge bookends
            if higher[0] - lower[1] == 1:
                # update existing entry
                merged[-1] = [lower[0], higher[1], lower[2] + higher[2]]
            elif higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                # update existing entry
                merged[-1] = [lower[0], upper_bound, lower[2] + higher[2]]
            # non-overlapping
            else:
                merged.append(higher)
    return merged
