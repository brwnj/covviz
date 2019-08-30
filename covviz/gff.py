import os
import re
from itertools import groupby

from .utils import gzopen, merge_intervals

try:
    from itertools import ifilterfalse as filterfalse
except ImportError:
    from itertools import filterfalse


def parse_gff(
    path,
    traces,
    exclude,
    ftype="gene",
    regex="Name=",
    y_offset=-0.15,
    track_color="#444",
):
    """
    Grabs the gene name from the attrs field where 'Name=<symbol>;' is present.

    returns:
        dict of lists
    """
    trace_name = os.path.basename(path)
    gene_search = list()
    with gzopen(path) as fh:
        cleaned = filterfalse(lambda i: i[0] == "#", fh)
        name_re = re.compile(r"%s([^;]*)" % regex)
        for chr, entries in groupby(
            cleaned, key=lambda i: i.partition("\t")[0].lstrip("chr")
        ):
            # apply exclusions
            if exclude.findall(chr):
                continue
            if chr not in traces:
                continue
            genes = list()
            for line in entries:
                if line.startswith("#"):
                    continue
                toks = line.strip().split("\t")
                if toks[2] != ftype:
                    continue
                # not currently converting 0- and 1-based
                start = int(toks[3])
                end = int(toks[4])
                try:
                    name = name_re.findall(toks[8])[0]
                except IndexError:
                    name = ""
                genes.append([start, end, [name]])
                gene_search.append(dict(n=name, v=[chr, start, end]))
            if genes:
                # merge overlapping genes
                merged_genes = merge_intervals(genes)
                # update gene lists to semi-colon delimited string
                for interval in merged_genes:
                    interval[2] = ";".join(set(interval[2]))
                # make the trace structure for the gene annotation track
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
                    name=trace_name,
                    tracktype="gff",
                    connectgaps=False,
                    # hoverinfo="text",
                    showlegend=False,
                    line={"width": 10, "color": track_color},
                )
                traces[chr].append(gene_trace)
    traces["genes"] = gene_search
    return traces
