import os
import re
from itertools import groupby

from .utils import gzopen

try:
    from itertools import ifilterfalse as filterfalse
except ImportError:
    from itertools import filterfalse


def parse_vcf(path, traces, exclude, regex=None):
    """
    parse a VCFv4.1 file, placing squares per variant
    """
    trace_name = os.path.basename(path).partition(".vcf")[0]
    with gzopen(path) as fh:
        cleaned = filterfalse(lambda i: i[0] == "#", fh)

        info_re = None
        if regex:
            info_re = re.compile(r"%s([^;]*)" % regex)

        for chrom, entries in groupby(
            cleaned, key=lambda i: i.partition("\t")[0].lstrip("chr")
        ):
            # apply exclusions
            if exclude.findall(chrom):
                continue
            if chrom not in traces:
                continue

            if not "annotations" in traces[chrom]:
                traces[chrom]["annotations"] = {"vcf": []}
            if not "vcf" in traces[chrom]["annotations"]:
                traces[chrom]["annotations"]["vcf"] = list()

            x_vals = list()
            text = list()

            for line in entries:
                if line.startswith("#"):
                    continue

                toks = line.strip().split("\t")

                # not currently converting 0- and 1-based
                x_vals.append(int(toks[1]))

                # info = toks[7].replace(";", "<br>")
                info = toks[7]
                if info_re:
                    try:
                        info = info_re.findall(toks[7])[0]
                    except IndexError:
                        pass
                info = toks[2] + ";" + info
                text.append(info)

            traces[chrom]["annotations"]["vcf"].append(
                [trace_name, {"x": x_vals, "text": text}]
            )
    return traces
