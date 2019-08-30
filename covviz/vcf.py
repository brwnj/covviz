import os
import re

from itertools import groupby

from .utils import gzopen

try:
    from itertools import ifilterfalse as filterfalse
except ImportError:
    from itertools import filterfalse


def parse_vcf(path, traces, exclude, regex=None, y_offset=-0.15, track_color="#444"):
    """
    parse a VCFv4.1 file, placing squares per variant
    """

    trace_name = os.path.basename(path)

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

            trace_x = list()
            trace_y = list()
            trace_text = list()

            for line in entries:
                if line.startswith("#"):
                    continue

                toks = line.strip().split("\t")

                # not currently converting 0- and 1-based
                x = int(toks[1])

                info = toks[7].replace(";", "<br>")
                if info_re:
                    try:
                        info = info_re.findall(toks[7])[0]
                    except IndexError:
                        info = ""

                trace_x.append(x)
                trace_y.append(y_offset)
                trace_text.append(info)

            trace = dict(
                x=trace_x,
                y=trace_y,
                mode="markers",
                type="scattergl",
                name=trace_name,
                text=trace_text,
                marker=dict(size=10, symbol="square", color=track_color),
                tracktype="vcf",
            )
            traces[chrom].append(trace)
    return traces
