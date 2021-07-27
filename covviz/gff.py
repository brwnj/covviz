import os
import re
from itertools import groupby

from .utils import gzopen

try:
    from itertools import ifilterfalse as filterfalse
except ImportError:
    from itertools import filterfalse


def parse_gff(path, traces, exclude, ftype="gene", regex="Name="):
    """
    Grabs the gene name from the attrs field where 'Name=<symbol>;' is present.

    returns:
        dict of lists
    """
    trace_name = os.path.basename(path).partition(".gff")[0].partition(".gtf")[0]
    with gzopen(path) as fh:
        cleaned = filterfalse(lambda i: i[0] == "#", fh)

        if regex != "Name=":
            name_re = re.compile(r"(?:%s([^;]*))|(?:gene_name=([^;]*))|(?:Name=([^;]*))|(?:name=([^;]*))" % regex)
        else:
            name_re = re.compile(r"(?:gene_name=([^;]*))|(?:Name=([^;]*))|(?:name=([^;]*))")

        for chrom, entries in groupby(
            cleaned, key=lambda i: i.partition("\t")[0].lstrip("chr")
        ):
            # apply exclusions
            if exclude.findall(chrom):
                continue
            # don't include genes whose chrom is not being plotted
            if chrom not in traces:
                continue

            genes = list()

            if not "annotations" in traces[chrom]:
                traces[chrom]["annotations"] = {"gff": []}
            if not "gff" in traces[chrom]["annotations"]:
                traces[chrom]["annotations"]["gff"] = list()

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
                    name = [i for i in name_re.findall(toks[8])[0] if i][0]
                    name = name.strip('"').strip("'")
                except IndexError:
                    # just grab the first item in the semi-colon delimited list
                    name = re.findall(r"([^;]*)", toks[8])[0]
                genes.append([start, end, name])
            traces[chrom]["annotations"]["gff"].append([trace_name, genes])
    return traces
