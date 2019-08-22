import re
from itertools import groupby

from .utils import gzopen

try:
    from itertools import ifilterfalse as filterfalse
except ImportError:
    from itertools import filterfalse


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


def parse_gff(path, traces, exclude):
    """
    Grabs the gene name from the attrs field where 'Name=<symbol>;' is present.

    returns:
        dict of lists
    """
    include = traces.keys()
    gene_search = list()
    with gzopen(path) as fh:
        cleaned = filterfalse(lambda i: i[0] == "#", fh)
        name_re = re.compile(r"Name=([^;]*)")
        for chr, entries in groupby(
            cleaned, key=lambda i: i.partition("\t")[0].strip("chr")
        ):
            # apply exclusions
            if exclude.findall(chr):
                continue
            genes = list()
            for line in entries:
                if line.startswith("#"):
                    continue
                toks = line.strip().split("\t")
                if toks[2] != "gene":
                    continue
                start = int(toks[3])
                end = int(toks[4])
                try:
                    name = name_re.findall(toks[8])[0]
                except IndexError:
                    name = ""
                genes.append([int(toks[3]), int(toks[4]), [name]])
                gene_search.append(dict(n=name, v=[chr, start, end]))
            if genes:
                # merge overlapping genes
                merged_genes = merge_intervals(genes)
                # update gene lists to semi-colon delimited string
                for interval in merged_genes:
                    interval[2] = ";".join(set(interval[2]))
                # make the trace structure for the gene annotation track
                y_offset = -0.25
                x_values = list()
                y_values = list()
                text_values = list()
                for interval in merged_genes:
                    # start
                    x_values.append(interval[0])
                    # end
                    x_values.append(interval[1])
                    # gap
                    x_values.append("")
                    y_values.append(y_offset)
                    y_values.append(y_offset)
                    y_values.append("")
                    text_values.append(interval[2])
                    text_values.append(interval[2])
                    text_values.append("")
                gene_trace = dict(
                    x=x_values,
                    y=y_values,
                    text=text_values,
                    type="scattergl",
                    name="genes",
                    connectgaps=False,
                    hoverinfo="text",
                    showlegend=False,
                    line={"width": 10, "color": "#444"},
                )
                traces[chr].append(gene_trace)
    traces["genes"] = gene_search
    return traces
