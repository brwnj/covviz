#!/usr/bin/env python
# coding=utf-8
from __future__ import print_function
import csv
import gzip
import json
import logging
import re
import sys
from collections import defaultdict
from itertools import groupby

import numpy as np

try:
    from itertools import ifilterfalse as filterfalse
except ImportError:
    from itertools import filterfalse

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
COV_COLOR = "rgba(108,117,125,0.2)"
gzopen = lambda f: gzip.open(f, "rt") if f.endswith(".gz") else open(f)


def validate_samples(samples, groups):
    # unpack and flatten
    group_samples = set(
        [i for sublist in [v for k, v in groups.items()] for i in sublist]
    )
    samples = set(samples)
    only_in_metadata = group_samples - samples
    only_in_bed = samples - group_samples
    valid = True
    if only_in_metadata:
        for sample in only_in_metadata:
            logging.warning("%s is present in metadata, not in bed" % sample)
            valid = False
    if only_in_bed:
        for sample in only_in_bed:
            logging.warning("%s is present in bed, not in metadata" % sample)
            valid = False
    return valid


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


def pairwise(iterable):
    it = iter(iterable)
    a = next(it, None)

    for b in it:
        yield (a, b)
        a = b


def clean_regions(xs, ys, threshold):
    clean_xs = []
    clean_ys = []

    s_xs = set(xs)
    s_xs = [i for i in s_xs if isinstance(i, int)]
    s_xs = sorted(s_xs)

    for i, (c, n) in enumerate(pairwise(s_xs)):
        clean_xs.append(c)
        clean_ys.append(ys[xs.index(c)])
        if n and n - c > threshold:
            clean_xs.append("")
            clean_ys.append("")
    return clean_xs, clean_ys


def get_traces(data, samples, outliers, distance_threshold, slop):
    """
    identify which sample lines need to be plotted and join up the consecutive stretches

    data - defaultdict of row data keyed by 'x' and sample IDs. value is list.
    """
    traces = defaultdict(lambda: defaultdict(list))
    for sample in samples:
        for idx, consecutive_points in groupby(
            enumerate(outliers[sample]), lambda x: x[0] - x[1]["index"]
        ):
            index_values = []
            x_values = []
            y_values = []
            for i, point in consecutive_points:
                index_values.append(point["index"])
                x_values.append(point["x"])
                y_values.append(point["y"])

            if (x_values[-1] - x_values[0]) > distance_threshold:
                extension_length = slop
                distance_idx = 1
                while extension_length > 0:
                    if (index_values[0] - distance_idx) < 0:
                        break
                    try:
                        traces[sample]["x"].insert(0, data["x"][index_values[0] - distance_idx])
                        traces[sample]["y"].insert(0, data[sample][index_values[0] - distance_idx])
                    except IndexError:
                        # x_values[0] is the first data point
                        break
                    extension_length -= (data["x"][index_values[0] - distance_idx + 1] - data["x"][index_values[0] - distance_idx])
                    distance_idx += 1

                traces[sample]["x"].extend(x_values)
                traces[sample]["y"].extend(y_values)

                # append slop
                extension_length = slop
                distance_idx = 1
                while extension_length > 0:
                    try:
                        traces[sample]["x"].append(data["x"][index_values[-1] + distance_idx])
                        traces[sample]["y"].append(data[sample][index_values[-1] + distance_idx])
                    except IndexError:
                        break
                    extension_length -= (data["x"][index_values[-1] + distance_idx] - data["x"][index_values[-1] + distance_idx - 1])
                    distance_idx += 1

    # fix overlapping regions after adding slop
    for sample in samples:
        if traces[sample]["x"]:
            x_values, y_values = clean_regions(
                traces[sample]["x"], traces[sample]["y"], distance_threshold
            )
            traces[sample]["x"] = x_values
            traces[sample]["y"] = y_values
    return traces


def parse_sex_groups(filename, sample_col, sex_col):
    groups = defaultdict(list)
    with open(filename) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if sample_col not in row:
                logging.warning(
                    "sample column [%s] was not found in the header of %s"
                    % (sample_col, filename)
                )
                break
            if sex_col not in row:
                logging.warning(
                    "sex column [%s] was not found in the header of %s"
                    % (sex_col, filename)
                )
                break
            groups[row[sex_col]].append(row[sample_col])
    return groups


def identify_outliers(a, threshold=3.5):
    med = np.median(a)
    mad = np.median([np.abs(i - med) for i in a])
    # https://www.ibm.com/support/knowledgecenter/en/SSEP7J_11.1.0/com.ibm.swg.ba.cognos.ug_ca_dshb.doc/modified_z.html
    if mad == 0:
        meanAD = np.mean(np.abs(a - np.mean(a)))
        divisor = 1.253314 * meanAD
        modified_z_scores = [(i - med) / divisor for i in a]
    else:
        divisor = 1.4826 * mad
        modified_z_scores = [(i - med) / divisor for i in a]
    return np.where(np.abs(modified_z_scores) > threshold)


def parse_bed(
    path,
    exclude,
    ped,
    sample_col="sample_id",
    sex_col="sex",
    z_threshold=3.5,
    distance_threshold=150000,
    slop=500000,
):
    bed_traces = dict()
    exclusions = re.compile(exclude)
    # chromosomes, in order of appearance
    chroms = list()
    samples = list()
    plotly_colors = [
        "#1f77b4",  # muted blue
        "#ff7f0e",  # safety orange
        "#2ca02c",  # cooked asparagus green
        "#d62728",  # brick red
        "#9467bd",  # muted purple
        "#8c564b",  # chestnut brown
        "#e377c2",  # raspberry yogurt pink
        "#7f7f7f",  # middle gray
        "#bcbd22",  # curry yellow-green
        "#17becf",  # blue-teal
    ]
    sex_chroms = ["chrX", "chrY"]
    sex_chroms = [i.strip("chr") for i in sex_chroms]

    groups = None
    if ped:
        groups = parse_sex_groups(ped, sample_col, sex_col)

    with gzopen(path) as fh:
        header = fh.readline().strip().split("\t")
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter="\t")
        for chr, entries in groupby(reader, key=lambda i: i[header[0]]):
            # apply exclusions
            if exclusions.findall(chr):
                continue

            data = defaultdict(list)
            bounds = dict(upper=[], lower=[])
            outliers = defaultdict(list)
            chrom = chr.strip("chr")
            chroms.append(chrom)

            # capture plot area and outlier traces
            sample_groups = None
            if chrom in sex_chroms and groups:
                sample_groups = groups

            for x_index, row in enumerate(entries):
                if not samples:
                    samples = [i for i in sorted(row.keys()) if i not in header[:3]]
                    if groups:
                        valid = validate_samples(samples, groups)
                        if not valid:
                            logging.critical(
                                "sample ID mismatches exist between ped and bed"
                            )
                            sys.exit(1)
                if sample_groups is None:
                    sample_groups = {"gid": samples}

                x_value = int(row[header[1]])
                data["x"].append(x_value)

                for group_index, (gid, samples_of_group) in enumerate(
                    sample_groups.items()
                ):
                    if len(bounds["upper"]) == group_index:
                        bounds["upper"].append([])
                        bounds["lower"].append([])
                    sample_values = []
                    for sample in samples_of_group:
                        v = float(row[sample])
                        if v > 3:
                            v = 3
                        data[sample].append(v)
                        sample_values.append(v)
                    # skip running test if everything is the same
                    if len(set(sample_values)) == 1:
                        bounds["upper"][group_index].append(sample_values[0])
                        bounds["lower"][group_index].append(sample_values[0])
                    else:
                        # indexes of passing values
                        passing = identify_outliers(sample_values, z_threshold)[0]
                        # remove those indexes from the list
                        for j in sorted(passing, reverse=True):
                            sample_values.pop(j)
                        # from remaining, grab upper and lower bounds
                        upper = max(sample_values)
                        lower = min(sample_values)
                        bounds["upper"][group_index].append(upper)
                        bounds["lower"][group_index].append(lower)
                        required_deviation_from_bounds = 0.3
                        upper += required_deviation_from_bounds
                        lower -= required_deviation_from_bounds
                        # trace data of outliers
                        for sample in [samples_of_group[j] for j in passing]:
                            # ensure that this point falls at least slightly outside of normal range
                            if data[sample][-1] > upper or data[sample][-1] < lower:
                                outliers[sample].append(
                                    dict(index=x_index, x=x_value, y=data[sample][-1])
                                )

            # update the outlier traces
            traces = get_traces(data, samples, outliers, distance_threshold, slop)

            json_output = []
            # add the area traces
            for trace_index in range(len(bounds["upper"])):
                for bound in ["lower", "upper"]:
                    trace = dict(
                        x=data["x"],
                        y=bounds[bound][trace_index],
                        fill="none",
                        type="scatter",
                        mode="lines",
                        hoverinfo="none",
                        marker={"color": "rgba(108,117,125,0.1)"},
                    )
                    if bound == "upper":
                        trace["fill"] = "tonexty"
                        trace["fillcolor"] = "rgba(108,117,125,0.3)"
                    json_output.append(trace)
            # add the sample traces for the outlier plots atop area traces
            for trace_index, (sample, trace_data) in enumerate(traces.items()):
                trace = dict(
                    x=trace_data["x"],
                    y=trace_data["y"],
                    text=sample,
                    connectgaps=False,
                    hoverinfo="text",
                    mode="lines",
                    name="significant",
                    # include color as primary colors occupied by area traces
                    marker={
                        "width": 1,
                        "color": plotly_colors[trace_index % len(plotly_colors)],
                    },
                )
                json_output.append(trace)

            bed_traces[chrom] = json_output
            logging.info(
                "highlighted points on %s: %d"
                % (chr, sum([len(j["x"]) for i, j in traces.items()]))
            )

    bed_traces["chromosomes"] = chroms
    return bed_traces, samples


def parse_gff(path, traces):
    """
    Grabs the gene name from the attrs field where 'Name=<symbol>;' is present.

    returns:
        dict of lists
    """
    include = traces.keys()
    with gzopen(path) as fh:
        cleaned = filterfalse(lambda i: i[0] == "#", fh)
        name_re = re.compile(r"Name=([^;]*)")
        for chr, entries in groupby(
            cleaned, key=lambda i: i.partition("\t")[0].strip("chr")
        ):
            if not chr in include:
                continue
            genes = list()
            for line in entries:
                if line.startswith("#"):
                    continue
                toks = line.strip().split("\t")
                if toks[2] != "gene":
                    continue
                genes.append(
                    [int(toks[3]), int(toks[4]), [name_re.findall(toks[8])[0]]]
                )
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
    return traces


def parse_roc(path, traces, samples):
    chroms = list(traces.keys())
    chroms.pop(chroms.index("chromosomes"))

    traces["roc"] = dict()
    data = defaultdict(lambda: defaultdict(list))

    with gzopen(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # header row is repeated at chromosome breaks
            if row["#chrom"] == "#chrom":
                continue
            chr = row["#chrom"].strip("chr")
            if chr not in traces.keys():
                continue
            data[chr]["x"].append(row["cov"])
            for sample in samples:
                data[chr][sample].append(row[sample])

    for chr in chroms:
        traces["roc"][chr] = list()
        for sample in samples:
            trace = dict(
                x=data[chr]["x"],
                y=data[chr][sample][1:],
                hoverinfo="text",
                mode="lines",
                text=sample,
                marker={"color": COV_COLOR},
            )
            traces["roc"][chr].append(trace)
    return traces


def parse_ped(path, traces, sample_col):
    table_data = list()
    with gzopen(path) as fh:
        header = fh.readline().strip().split("\t")
        datatable_cols = [dict(title=i, data=i.replace(".", "\\.")) for i in header]
        table_data.append(datatable_cols)
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            table_data.append(row)
    traces["ped"] = table_data
    traces["sample_column"] = sample_col
    return traces


TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8" />
    <meta name="author" content="Joe Brown - Base2 Genomics" />
    <title>covviz - Base2 Genomics</title>
    <script type="text/javascript" src="https://code.jquery.com/jquery-3.3.1.js"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/v/bs4/dt-1.10.18/sl-1.3.0/datatables.min.js"></script>
    <script type="text/javascript" src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js"></script>
    <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/v/bs4/dt-1.10.18/sl-1.3.0/datatables.min.css"/>
    <link rel="stylesheet" type="text/css" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css">
    <style type="text/css">
    .disabled_div {pointer-events: none; opacity: 0.4}
    .selected {background-color: rgba(161,234,247,0.5) !important; color: black !important}
    .tab-content > .tab-pane:not(.active),
    .pill-content > .pill-pane:not(.active) {display: block; height: 0; overflow-y: hidden;}
    </style>
</head>
<body>
    <nav class="navbar navbar-dark navbar-expand-md p-0" style="background-color: #3C444C;">
        <a class="navbar-brand m-0 p-1 text-light"><img class="pr-4" src='data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAASAAAABzCAMAAAAPMC9fAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAHCUExURTtETZXt+0NLUwC8xP7+/gC2vjpDTP///wC3vztFTkVNVUFJUjtDTDtGTzlIUTtHUAAADAC5wT1FTj9IUJPr+Dw+R0AvOwAAFDlHUACstAABIkAtOEZNVQCepgAeLz1HUDZASQC7w4mMjz4zPQDJ+wASKgC9xT40Pz07RD44QjM9R8nKywCwuNra2wAZLT5LVDxASQCor/v8/AAKJ+3t7QCZoACzuwEAHI6Qk5SXmerq6pjz/ylpcvT09DlLVNfY2CQyPgCjqsrMzQAAASJyei1iajFbYyo2QdDR0sHCwz42QC86RPDw8fj4+ElQWJbv/DRVXR4uOwDI+TZQWACSmhp4fw8mNZudn7S1twCNlOHi4iwDGA9+haqsrgDO/7u9vgcjM3N2ehcqOE5VXJfx/gCIjwC+xt7e34CDhsbHyFNaYeXm5qGjpWJnbGpuc9PU1TQuOFtgZo3f7Hp9gTY2QJTs+YfW4jAeKwCFizElMc7P0AWBiD0sJ4DJ1UEkMZHm816Mlm6qtC4TI0pncIWNkHS1vz8gGFeBimSZonm+yURZYVJ4gADByRykygDC8U9weAC65gC4wCSWuS2Cnw+x2zlHTx683ycAAB05SURBVHja5JmNT+LYGsahl3uatiikxWWDXdrokLBMpAZMUDoqkdnpOCAMs7PItUEFyfrtIAIEI1E0fmfGycze//eew1cpIA7ejEnwjUmkpe05v77vc573oNH0EgZBEDhgz6aLBRTFdJYfQIcIzfMPBIfIFuTTNHdzd3lSicvbM7JyiIeQNPrnTYfLnsr2s5PydT5zbNKGUykxlQrnjjP50vnVrVCQ0xqYSM+UESEIWbl4dnWdP06JYiqszeVM1chBUqIo5jKl8u24XCQEwzPEo+eIQuGsXIJwUtoal0bUPodFUZs5v+TlLCy1Z8cna78qmapwcvdEFVI4Xz47TQvM8wLEZW/y3ekokMLicemyUBS456RFXDovtsExmVprrXZYK2ohorQwoH8+fC7DplY0UJvD4VQ4rNXm2jVJDJduZV54LoCE06uUSZm9FsoxrLfjTCafz2Qyx7mUWJPuJkZi7vrmVDA8jyTiiidVQBWJSZnypfLV5d2NQAwNjQs3Z7fIGB2HK5AURDkxc3U6JPQzIWJgoPaPBmREU4WONnN9dTZ0Kp8Ws0M8IOAZPpsuyHL65rKcPxZTOYUR1KLSTUHo3w5k0GAfH6zWCJe+PYZeMAVtDpBPswRqvQwDEA9BEAMGDn3ki8hGQiugpBGss8yJzPVrmc0cvXvz55G9TuimXLo+0aBmQjB0SooBDjUixbvzjBhWEIXD5QLoywWfMFx8//rly+c37yqEAGeXT+UC6N6zo2ZNJssZmEVKEpXGswLoO+c8YDj6trz84sXy1/dDVSFCOxrCQ42oXsMJvExARLkGIjF/k+43QnoDcfR5+QOK5c8XM7WD+h+qFD0h6GXyPKf4AjFzVuwvQpDPMOTzAsaH5a92otd1iBCG5Nt8I4lMqeO7/iIE66vGBwEa6hkQrFAhnS5r6+7bFO4vQsSgwqepxHpDxBjkS+id6oQyZ+m+sYz6mYtvNT4wgb78YX/c/hchFG7yDUKpDGnvk9Ue8vnexOfV6OBj5yVk09cNQmI+qxnoDz5HryGZugC9GZ559HvXcxr5vL5LYhKv5e7NPUGrAyrf08sW0NAMhcGgGPoer2d//0UR6D9HFT4M1jWodgupH+DkspJDV3JPMkRQJEs9cSPHsBiQoja325mUNCTWiZHh6B+Fz8zvdT6AwSSnrUs4owGaxIg2IWoQMmlNt+luOcQ7VMEHoh5zlCafEBGFRcxuam9janFx8mCX91gcJAPuFWioP+/fNfhgkmVvcqpLTG7sSkYbIFuhC6jKakKd5wFx7+Csc0GvOtbfxjYkD009UaER5Ko7EZudWMBxHfxb8L7d2TIGSPXTDaOvPlQF6MOLL6+GFT7R3VkX3j0WvC+naCNPtvawhVKdkHh+vwxhlvUO99R541ISewpCgCHMG2Mtc5ze3Ldgqq/9d7RRYC++HzX0h4okpuFou0bljsFYyIJxqglxYDxfbzu0d9n7fhDC3GO4bmzuUIlP60F0y+2EjQRPwMfh26y85vX5+OLk5GJ8c30BfdzxUKoEet1YwL5dNPgQtLSNP5RAeJVScMfCq4pCL6TPjqvdvUks36vTCBB+4HebG+FZDY3Mwae6Rpw/P4doEIAvCPfG961Gs9vptJiNq7txL5zWSnMODR59r7dg/xwp/hCzTeE6/EcCIfpE+zA1Ibm2o/0goJFICMPISmAkRvFW/94SfIu7vp/+KxtmmYNTjEU8AYqsrsokJZmlFTiluEfJ4HoGQQF602SgSeP8DwKqIPLuJ1sJVQyjyZQ6KQqaboCsquSjqVBkFR6edWM/m49zEs5w0U+RTRLAsYx/EWbwXtP74ZAGoU2OJgHSANL48ocBwRviwX2faqIDvFDZ0hZLaQPxAKAWaWADIah+G1HqJwNCT1/xr7WObs0Pc+ilUVl4DL//8XUZxrfhwaZR9gYIEvLSPKey1NmzPOztS8KQQdMTIA0I+WM4Pm9kfyofLrC/gE9vObh2aSKDeJDlm3aC3hm+f/72+kjVNvUICBJ6a1EVGeCy6duTuwJ///79PYA0jC/hwsdsVEtDwDGoFwD3dAsUhVEUQ3Q6DVArQVHqaynrCI6PdSpk1gy1aSPavFemP7q4GB1U3bNXQKic1f5BzxnSxSzX5bfo+wDRfCiIb0u0kvyAwEg6EElGHAxb/TrNcE3rNUk5Ilbo7FGv0MIIQLPPw5PWiINqNsmYBWrNpplsHxZpjuNwMurdoJkZdSE8BlCQaTXNHMM9qAIdAMEpTuATQKPcDaOdxsB+4iCxFTFGMSR2UkRS7LBktG4lRpCz5802Qp0UJO9xr8GTI4k1q1EiG4aNcGwdHOw6OrRepGcHx+PmB0q0Z0Dw2/FO7+MRgGh+bRrfjjBEvVgxtzQ5tz3t0rmml+Y33BGW8c1NfOIrBJEd3o2tT7uqzv5wymHBaNDILcbM7swG0UlXcD22Z3ZgCqFkMtCpNSU9MIN2LA+4zMcA8jZXxf8BiEpu4PhhnTageMuit+k5YwlPKLmNB3kNXbHD0spC8zi8i5ZALU8AFliNB5tbo801Z8Ph0EiWOgRr3ETDekQGtbcZLYQ6TLZ3QIBe889BEajpJ8B8zCeUAEubsXis0gy44m7fEu6tAKJ5DWrovJvxKdgsxA4RjsMtnq5emtyqntwZGZncmV+Cg54esTzQxcBLg/j0mqO2xg8OGn4wg1qZ6NoB9bg0VwGtYoQSNEOxNLJqdY0GlLQG08c1n0gazR7UDMSmkX2pA8Ish5DAVMDssVlhs2AmdyZgfdiwyqWRXchrYhEYLVar1WJMbszCYU51JQToEPJBc9X8JYjR4eHRjr9gtAFC9Jdm67E+0Y5Ihy89JoOcFFCCl6xGENOhZba2WgGwDQkk/D6KJUmSpST/PrxqaqwKiIrCatwGFoYlMQq1LMBMzK5YKtdyDhaO8pC3ECxGURTGUlFjXIe7EpHOowQAEAwL/DsNJ02M2//6+PEvu554GBCc/dRa0lIPJzh42UoIrmMhnn5EBjXtua7tTa7At+6adNfElPTAcltyREMMnAD6Y0KS+y0e9OITEBBgUUd0YAvVLQ4g2MCqD9QWcliaL/2ArZ0EGopFybkUZVo3fitbriRLcj4jMw+nUu3FCD3x279g/KbvQKgFECwfj41XdmApTdI/stCaY649iekRkC440RyV9WZhNmGr8cGQnQtiUlNVAAwkl+DDEKAKBC+gVQaVqi5/qN3Gl5yMaoMw5J9X5K1KB2aWw7dqc7s9RrMvEUMiFqs6ul9/+fj3v2H8/Z9ffn1Ag3T4ph8jOZIla8FSVMg/0qbTB0mqV0BtSubSzfF+tjYtzAaTbNLCqjw6Fk38j5pr7Ulj3cKTKb44DAMjBEanIyDEFFAJCnKLbFC0VkFSogIlJQBmbxRUYvzS2OT8B//xWe9cYK6oOefD7kqatg5z4XnX5XnWekebCqA0Ml0VugkwPutab8E4MAgVf6YcVIJtYp7Uat23XyZ4fYotTjxOOf/afl1aXl563XYZXEgLEDxNlfB6UJZVrFL2opuKQfB/TGFKHlRIpwvSn4LsQfaLXrkqSl+6AaJjlNJlVSRUIDPLIYZLcuXGa4IPZKdvnKA/FU4oPlWDMoPwZOO5sw7cdfZFimdjXiIYTOwA47O8/Jr57qMWhhjmgIIn7un2ZjYY0sHqgz7IPu5B9vCmJzQ3AXIQzv+dm1NACDm4NuZsesEkBp7kQYkchOAjrBck6SDOUXMkoBjlIVQRocyScBITAy+XICVXbJS7abWgTE/vx2xcXg+X//cnDNDSytpbHiRmdW9zqqZckG7oQEebpy4ePpyDdDwodJpg4/dQyAteyPcIR5DJRSlE1iSAQHNOsccNymy2XwXFRissmk7sgKsgh25SJVTBJyXKj8jA40R0ms5Lr9vttsKPwWw2LihqxOXflwDaWI+4F+YgG+ZNJLjsnCbuAedxkNyLFiBQY57/Aw9in4D5vGRhlXHftxA38hCcmiQeFAwRmNzsTe66tw9xng2EpPQF9CBt71QC+lFVgPMW7T0WnATRzQGm4KXWU5MDjpXlA2VEzmUKwUQPN0SAPv347FpYxQCgsafezGuyNlApXdfRZt49+DCTJuoBiF37bZ+mEJ22j0zIlYM7kgECqdG/n0hPUZxM22ESFKmYjcdFe7G0UzLYzgVeXAQc3AEVqxbOpqo07vfiyaqmE+B0ra1IHnQdWX3Lg06zUNY1AG2ZANRmHf87QErHDFZA7HvsbJImy3cnAwRl3cs3btultDxhqN094qGE2A+zsrMK5G789PbamBdIC6/3fc28Lot2EmXe4EF7o53RzkTDeYwAfbTKWwJEiOk/TaEgCoIHNc086EUBCF/H22R5z2Ouu9Up4oftbtJSPa+pR0pzK+X7NOHBH7CH+bql7oAihtFZMi1ib6l5E4DUSoNyr+qM+ghAhLc8kQpDo4OBoowz2ZEKIIKiBRKVExw3fGrtSESPIkIFyEF81sS4MuRhsZyNeMFal305XpGL2Oo7mLROwO/ZjyoOB6clk7OWojsZ8WstknR/BCBcgoAzeLGn7A0MEyBx9dNzgKQdGqAVvNUU14JYA7mF8/jFA5pRW5WJaVxgIYzbizpY/p9Sjl7Z9b1DixmEfT5AkqmSmit1hvIXcZ//df1Ta9eXEeYjAPVHWFR6Sdy7gqKM9EUsb58BpGrMIiJI1rHYhCws8iBYMmRiUtfnzIxiqQH6JVf5wxjzoYYZhqeTa4LiATYyT+S47Mj4HG9vbGxgmfdp4xP+G/6X+XnufjdAFEHXYPmrXu8pJKNC2TAYypZmAOk39XgISOwlngQmbYMQ0mOLaHk6JXLwxQD9bU2DFngQ/uEoFwjQqJ5tqR2ox0rxzERPXleW9LaycWxAyLqK8eAhkzIsAY4xe7eiyaRwOKeIVcI7bBAaOUoRVEGsfHS/A36e1SZh5Gg2SfQ+gD7/UADyu9/tQfhHpXAiQZPeeoJMqxxo2pcDbFUmoFoDR30vQMF6uQpMsYc5gyjGLgabqq+JhHK9NlPzlW4phdSC1BuH+j7NkkhUIcXxUH2qt17J5yqiD70foE/7RjFvDhCOLdv0lofs4xVO2cFEhc/FuCE/JuXbXTEBCNcCi9k8PTNxPl6t4EapNNNDAlZjtUcWb/3C2cNDCptkxz6ZYIAochNy1Ld+X5COIkTRIocSCZnUbBxnVadWK3mb/V5MBRJACdXNVUbSaoDMPIgwb7nuvTxnm8A4AZ7HI7umhJWG8hSCiR2aA0SZAnRbGQZmlkhxbEBsqmKdKRc0LLZaqVQc8106FGDD4LgDUWpQdOMW6ODkmUuIRx1kqF9p2ex7onzzhDzgiMU8n6pChSO98QT7tAVXbs0BalUaATNL9NEbOcjgQRiei7tHruEQ4Xm427PrOo5KDnL5f/2jST8yQLtO0yTdvc2FZ5Zrdc/EPTq2e2VOG0RlLLZG+XEgy2UD9dzUZrflKxORByEyixGyTfM3AY5luYQQnooNQXwyostjTG53Wk9DPsuXx7kj3PHPV8hZiPXUN1dZLgyfWVjFjB1Fe/HsgS07vKbwiADOsolv/Z8VlSkAGdgEuVmy2II0Gszm2MgbCvREpTV62fombrCq3W6eykQRORJ4twz8bHR0dna0I854evLUGpGn6EWS66VpaYL90t4ZyNwQl3nr/U8XT9XFPMjgQcX2mD0VAJ4y+yjCY0hQeNQnJSGn7/fh7szWXpct1oHkO2ZT/vRRODGcl2cQW9yzajPgxVk9UO8XlLkYWW7ea0RQJzzbIIboUDY3UiWCTqvZkK8s4NmStVIAgD7CpG1tgY0LQdrRYB+3zODRMGmK8J/7ZHP+xzqSvY3ut60jtd2dtVsDgg3Qmq00lGPID9o76WKxWJh2n9gQ6TltT3vSbBrRWKz2SoVasVgrlNrPgYBDNcQhE4nndmlSK9ZA6j/3+dmV6X5Ld3PNg9yEPqDFsJo/rXs8ZJMbvFjAI1LpvrJ0LkZp8My6Tp+u/C7DDCrAchqDNMI3kMOwIYYmG2wKjZ/GIZ6rOjxQkRKV1NwTQaym4jdPTzdVnh3S2i40SQ9Zvlq/EU6zHBxTNfeHFc7Ssgir+SWp5Wqm5vUAjYHFV/nnbzYreKzUvDumdJ3+NnadEKmTSoIgkObzYLzdO16NE6RDimNSUE8mQKwS8Wo1jkcURt0CJCAUD3lILfCINtNp8nM4xP3RC/pBhKGjSFcfSmbz1MX9IMrJLLqLhU4y3wJEeYIeSvmA7pPiUThsun0Iia88UIZrowUmSgKJ0L2ro3gTd6Tyi5xH6igaelsM+KlEhF4PYn/Yiy2qnvQX99sAYea+Z75Tet6TpvU9abdSCl4zl0nqTwPol/VUw9CTjgvgQebabNFUw6XwdfOmyr/amNiBLAnMvN+Yg8qDuzOd9TpvzcWYWY7e97v+MICo5KU4WYWH/2WswAaAKE+1wmpts9LTJW39ZJXyfc9I99j4+ccBRDDn8upuHJ4zbwLkQR59WbzZ1M+e9R60GrneWLKeDPzbk9BnieQuv2a+GhKEPgcFBeP7dDcp3djHsLtD6Q1BnvsDX15lZC4NMfbb4P86gMpcijdYXxtixv1Bqii++mwZYRTjZgxv5WEjKE2nUP+P2enMe9O/8S7GG7tXV+eu7lRI3KHe//Wz+ZKpTdNagPRD0FVZEIMHHX+x2HDvdvti51Gn660v6TR/lX+VSUajSbdbh67ydbDksXSQVacv6qNcGniS55FIlHArMXYlN0dXTvR17J27XA0tIfMsZ8ok5A7/l+8HByeuz8k3fvsZZbYRjmIiscuvXy+/aOD3+aLyglNRXzRpccHVqN/19eS7zz97M9vJJP2XB9fHJ0l5BKOKsX1dBCCzuZjVK3VWu1yZ2Ily/SuLGsZEDtYzyytLa1eXSqGgkuvba6Ltbv+S5QlzfpXZjwFAztXD7atZSaGcsevD7aXM9vpxRIkfEDfru7vX4v3ckR+7uz9k8JjYVeZwvgqM/2R/LbOU2f3tU8KH8V3ub/+3vKt7cd244iITXaRcadQhaCpVKM1EoA8QQYIrJYAQNUgKqBRcbNGa7lPAhhLaEvpisuDHPO2LgP6/PWNbI/lrc3cboHhngV175RnPb87X7xzNiCAmw8fVk/TxJaN/7X3SV2N1kKCivapherAmzCcK8VFeDKtu2ohsePOfBstleK6/Z7sajp9E3kTVNJchn7ENIstwBCiGOTahoWHaQU8uFX08jZkdw9v2zPeR77P06KLUpLR9Jsc5Y2h10KghmL7Cx36FO+0FD7sZRoONQr28rIoq6+3Bk8KqsThzM2hudxQWI3AJ2DGKNWz7ghnpQeYTe11Vq5j4q2NmndPjPucjGjwIy9kQoUIfvgDI8mDkfrZ+XKUoPo6M6dpXsmYx72zkaoOXya8nhV6zV+ODdFpEHwXoipccnBzJi10btk7VCLkw083MiVpowfCdYHIMvmMT6tYIEOd5ZIadINjNbTIEWvv8gaLAipcQxfdjCA99bIbKiprMcyJ3lFInWQ8jWzTz4xBEx2mqQDsz0+d++FfY7aOH9VDQiMurRBXTR0Yed6WeJBqNhqOK9gB5h3zkyHldlMss1aYSBDqD8iYqDaukhTJMgAPkz1JCqr+tGHHtzTWAQHUY2e4wtrBEB/fAFzSvAkrDUGy91EPQgmuh7ruX7xf76ZvTVB58f3ZDf4cr+JKVmqTOG2hzQU9SP60L/lMMWyD55FIwV5lTCoDAZMT+cfZGCJ85LAIHaJM1PbE7cJ7zC4DUQd3i8mDXxY0nVlshIs+W22oeCI85qsFZ5fPlOw7fnaY6MO2O+NzMdKghyEpoJAtbznNZiP9YeCQDi95PLnIZ23oCIFjc3F8fZg/T2NgHMT1I0M9AEkDRalVgeAKQlMDIFwxdj7YKAt9A4iUevvFehI61YaCTr9+z+uFszypfX/K8AB2UKQKAZJ8RESpxgJRDwQhNAPJtyQLjUqUbAVAtowlA8WIEyI0AzRxtd/MbAGnQzQUBldS2WbtpTMhmloyEbLBC3E7rEwl60a7nH0Tl+dLFgxapt9iyb2u6ZlVdkSEba6OKVYc2iN4eoLKtobN4VLHD7LVD5htUTHgxAEhbxE+ZJy2uA6TD9fYxchqZhYo1k3rJokoVBoHVGV06I/Uv2zf/05dn++aFgoEA3Ux0YGfF+sop9cj8ebuZAAReLIx400Yj7ds6BuaiKCIBqocpixdmqescusFbHQAKSlqlWNMuARIjk2pXWgY2PNMYXJszl3Cpq06tsLUzcSVomMpUF37P74L79BfB2Z+88JdvTk9egDC+yY+wk/RmnoPfKAPxD7j5cLf0TwHau/l24uZ5DMPFPRcFNos7QTeiUeToKSHHWPQIkK6aJlCHGwBxD9vH9S5oKX2o6VEiFyv5cReZya5S2NaRzm8H5pPpRjP0/s///LizOz789ZPJJshDl5omDG1fRzf3hRve0u/zVVV0rszs0c2zOMvOA8VDkAf2iYwp9CRFffpQFI9272eBMVExeMHzAWcATQJFuoKRl1314JK43psWUFmZycuiqR/inhThdB0HO93LdYvF6S/f//ivXzi849M//fD3b7/87bvPzrMIgYuuSeVl+sHMEPMJD/qvUA00pRqcJkBwPo83gmrAqxQh1PcIIVc6Tp0D9DQbzNQif8quUg01zBiMTJjPlO6whJZX9QyRXEbMX03ZL3a2w3RIPA8FQu+/+90/vvr2mfbVv3/8w2++lj47FR9JNUSWgJMM6ZlExp5tcrIaZ3NBVs2BrNrxlg5kNctnXD2wV+WZMJW6qa9tue+BcYbJyMVmQGz1I0BpvvQmfQiyqkrtMHIT6ENmoZjlPSN5+hCdeROx4MRejAh9/sl3X3z9TPvi+z9+frFPneOzGroDkXz+QD1V90Ke7iipOZYVLTy00Xrp+DhNUx0vBA9IJYgnG5NO0z4W/s/k78s+RKbkOPJw9IFmRO28qKrGPLvjFGj+MWYBVzZBCCz1+ysHl4nf7/n+msspC3ye82CTjJcZtVEyvZFaG9v0PWHfpqNhLWzbUD0p2l79nHRx8KpxHHkChm6YUXDe3T7ILnpFIDQP/of6gyHRpcDHX35MKeNKylXdpwU/MnNsGPozqVX1uZTrlZGNy+542OI8DNUHsEO1g19Z5NOwaWaj/Mxefibs/2+FY9QLIneO+rojxbG3mIl+mI0TXbqXZtClL7w9W4bhy9WMHyhcxBNJnIfW3eDDY5dMICSjtHHUF85OxZG5ImTEp4nwPZ3brqqeQAjUQ96GLdZfAk/i1ClSJvi0+L7OtQeEwP8ME1RQWtDwYyECb+sslj0ZAbbnEb63c/9Vla6ZIuZISFbT0PgIRdNxQvVtjMbPorQM8f09FwFi4E4ZpUBGclYHnva8GKlYD51ya6NefFBBWWJad/lkDUxrW6gZn6nidroTadjSryOKdZNG9SpGRJ6I3trT7vUpUTjC7igLABFj9rJQKQ0lA1vGvtTN8wi6bmCsJq0TNY8zZQqPjOLCMe722TWalXgPOZNHiGSCiJ11jeQ5tA3NZE+REjPyALRFsZrJCOzWeDlhWRlg6Y6batDGZWQKkczr1ba7eijqeYlBfMpFU3TbLOWF7Ak6XCfjjibWnT8/C5tBx23uZOLwN0EIYMv3qZo4lnuGEOnlyUVwFVOWJbXu/+G0ukGtdX4C0R4keN0f6jH94fXJvxVG3NpJsPQWGtac+Sqf2N5TnC7eluUe9W7hhfitPNuYM4fF2mYnJuZG2xupPCuC1tClt9M4RGrn5oj18m2QDtaJpOuGRsZbe1ovECzPgzAnZwdzfGaR968JQwoESoljWrr0BptlmLSdd8s03zuxXgDDrTVDTLHdba05kYTf7HPVgWqZHjWbapvNbPDuRy8m53bqrroaezSUsCG96aZaWDUDSiN10dRFVXVVxSNG03NoqBlYl95+29NSrAPDCNsg8II2Ck0N3rL0O9as/wKiMOxXSvdq/AAAAABJRU5ErkJggg==' alt="" height="58">
            covviz report
        </a>
    </nav>
    <div class="container-fluid w-90 pt-3">
        <div class="row">
            <div class="col-12">
                <div id="chrom_selector" class="btn-group btn-group-toggle d-flex justify-content-center flex-wrap" data-toggle="buttons">
                    <div class="input-group-prepend">
                        <span class="input-group-text">Chromosome:</span>
                    </div>
                </div>
            </div>
        </div>
        <ul class="nav nav-tabs pt-3" id="myTab" role="tablist">
            <li class="nav-item">
                <a class="nav-link active text-secondary" id="scaled-tab" data-toggle="tab" href="#scaled" role="tab" aria-controls="scaled" aria-selected="true">Scaled chromosome coverage</a>
            </li>
            <li class="nav-item">
                <a class="nav-link text-secondary" id="cov-tab" data-toggle="tab" href="#cov" role="tab" aria-controls="cov" aria-selected="false">Proportions covered</a>
            </li>
        </ul>
        <div class="tab-content border-bottom">
            <div class="tab-pane fade show active mb-3" id="scaled" role="tabpanel" aria-labelledby="scaled-tab">
                <div style="height:600px" id="scaled_plot_placeholder">
                    <div class="mt-3 d-flex justify-content-center align-items-center bg-light text-muted h-100">
                        <div class="d-flex flex-column">
                            <div>
                                Loading...
                            </div>
                        </div>
                    </div>
                </div>
                <div class="row pt-3 mb-3" id="scaled_plot" hidden></div>
            </div>
            <div class="tab-pane fade mb-3" id="cov" role="tabpanel" aria-labelledby="cov-tab">
                <div style="height:600px" id="cov_plot_placeholder">
                    <div class="mt-3 d-flex justify-content-center align-items-center bg-light text-muted h-100">
                        <div class="d-flex flex-column">
                            <div>
                                Loading...
                            </div>
                        </div>
                    </div>
                </div>
                <div class="row pt-3 mb-3" id="cov_plot" hidden></div>
            </div>
        </div>
        <div class="row pt-3 pb-3">
            <div class="col-12">
                <table id="ped_table" class="table table-hover pb-3 display nowrap" width="100%"></table>
            </div>
        </div>
    </div>
</body>

<script>
const data = [DATA]
const sample_column = data.sample_column
const cov_color = 'rgba(108,117,125,0.2)'
// const highlight_color = 'rgb(255,147,0)'
const cov_layout = {
    title: "",
    margin: {t: 10, b: 40},
    height: 600,
    xaxis: {title: "Scaled Coverage", showgrid: false, range: [0, 1.5]},
    yaxis: {title: "Proportion of Regions Covered", range: [0, 1.]},
    hovermode: "closest",
    showlegend: false,
}
const scaled_layout = {
    title: "",
    margin: {t: 10, b: 40},
    height: 600,
    xaxis: {title: "Position", rangeslider: {}, autorange: true, showgrid: false, showlines: false, zeroline: false},
    yaxis: {title: "Scaled Coverage", fixedrange: true, domain: [0, 3], showgrid: true, showticklabels: true, tickvals: [0, 1, 2, 3], zeroline: false},
    hovermode: "closest",
    showlegend: false,
}
var ped = false
var cov_traces = []
var scaled_traces = []

const load_buttons = (arr) => {
    var btngrp = document.getElementById("chrom_selector");
    for (var i = 0; i < arr.length; i++) {
        if (i == 0) {
            btngrp.insertAdjacentHTML('beforeend', '<label class="btn btn-secondary active"><input type="radio" name="options" data-name="' + arr[i] + '" autocomplete="off" checked> ' + arr[i] + ' </label>')
        } else {
            btngrp.insertAdjacentHTML('beforeend', '<label class="btn btn-secondary"><input type="radio" name="options" data-name="' + arr[i] + '" autocomplete="off"> ' + arr[i] + ' </label>')
        }
    }
    let chr = \$('#chrom_selector input:radio:checked').data('name')
    build_cov(chr)
    build_scaled(chr)
}

const build_cov = (chr) => {
    // hide the placeholder
    \$('#cov_plot_placeholder').prop('hidden', true)
    // show the plot
    \$('#cov_plot').prop('hidden', false)
    cov_layout.xaxis.range = [0, 1.5]
    cov_layout.yaxis.range = [0, 1.]
    cov_traces = data["roc"][chr]

    let cov_plot = document.getElementById("cov_plot")
    Plotly.react(cov_plot, cov_traces, cov_layout)
    cov_plot.removeAllListeners("plotly_click")
    cov_plot.on("plotly_click", handle_plot_click)
    cov_plot.on("plotly_doubleclick", handle_plot_doubleclick)
}

const build_scaled = (chr) => {
    // hide the placeholder
    \$('#scaled_plot_placeholder').prop('hidden', true)
    // show the plot
    \$('#scaled_plot').prop('hidden', false)
    scaled_layout.xaxis.autorange = true
    scaled_traces = data[chr]

    let scaled_plot = document.getElementById("scaled_plot")
    Plotly.react(scaled_plot, scaled_traces, scaled_layout)
    scaled_plot.removeAllListeners("plotly_click")
    scaled_plot.on("plotly_click", handle_plot_click)
    scaled_plot.on("plotly_doubleclick", handle_plot_doubleclick)
    \$("#scaled_plot").removeClass("disabled_div")
}

const build_table = () => {

    if ("ped" in data) {
        var ped_table = \$("#ped_table").DataTable({
            data: data.ped.slice(1),
            columns: data.ped[0],
            // scrollY: '600px',
            scrollX: true,
            scrollCollapse: true,
            paging: true,
            pagingType: "simple",
            info: true,
        })

        // track if we actually have a table rendered
        if ("sample_column" in data) {
            ped = true
            // register table clicks on sample_column
            ped_table.on('click', 'tr', function () {
                if ( \$(this).hasClass('selected') ) {
                    ped_table.\$('tr.selected').removeClass('selected')
                    reset_plots()
                }
                else {
                    ped_table.\$('tr.selected').removeClass('selected')
                    \$(this).addClass('selected')
                    let sample_id = ped_table.rows('.selected').data()[0][sample_column]
                    highlight_plot_traces(sample_id)
                }
            })
        }
    }
}

\$('#chrom_selector').on("change", () => {
    let chr = \$('#chrom_selector input:radio:checked').data('name')
    \$("#scaled_plot").addClass("disabled_div")
    if (ped) {
        let tables = \$('.dataTable').DataTable()
        let table = tables.table('#ped_table')
        // remove all filters and reset table
        // table.search('').columns().search('').draw()
        // de-select all rows
        table.\$('tr.selected').removeClass('selected')
    }
    build_cov(chr)
    build_scaled(chr)
})

const handle_plot_click = (click_data) => {
    if (click_data.points[0].data.name == 'genes') {
        let genes = click_data.points[0].text.split(";")
        for (var i = 0; i < genes.length; i++) {
            window.open("https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + genes[i], "_blank")
        }
    } else {
        let sample_id = click_data.points[0].data.text
        if (sample_id) {
            highlight_plot_traces(sample_id)
            if (ped) {
                let tables = \$('.dataTable').DataTable()
                let table = tables.table('#ped_table')
                // remove selection
                table.\$('tr.selected').removeClass('selected')
                // run the search
                table.search(sample_id).draw()
                // highlight the selected sample within the search results
                table.rows().every(function(row_index, table_loop, row_loop) {
                    if (sample_id == this.data()[sample_column]) {
                        table.rows([row_index]).nodes().to\$().addClass('selected')
                    }
                })
            }
        }
    }
}

const reset_plots = () => {
    scaled_layout.xaxis.autorange = true
    // de-select scaled plot
    Plotly.react("scaled_plot", scaled_traces, scaled_layout)
    cov_layout.xaxis.range = [0, 1.5]
    cov_layout.yaxis.range = [0, 1.]
    // de-select cov plot
    Plotly.react("cov_plot", cov_traces, cov_layout)
}

const handle_plot_doubleclick = () => {
    reset_plots()
    // de-select in table
    if ("ped" in data) {
        let tables = \$('.dataTable').DataTable()
        let table = tables.table('#ped_table')
        table.\$('tr.selected').removeClass('selected')
        table.search('').columns().search('').draw()
    }
}

const highlight_plot_traces = (sample_id) => {
    let s_traces = []
    let c_traces = []
    let highlight_color;
    for (var i = 0; i < scaled_traces.length; i++) {
        // let trace = scaled_traces[i]
        let trace = \$.extend(true, {}, scaled_traces[i])
        // limit to significant sample traces
        if (trace.name == "significant") {
            // de-prioritize; gray
            if (trace.text != sample_id) {
                trace.marker.color = cov_color
            }
            else {
                highlight_color = scaled_traces[i].marker.color
            }
        }
        s_traces.push(trace)
    }
    for (var i = 0; i < cov_traces.length; i++) {
        let trace = \$.extend(true, {}, cov_traces[i])
        if (trace.text != sample_id) {
            trace.marker.color = cov_color
        } else {
            trace.marker.color = highlight_color
        }
        c_traces.push(trace)
    }
    Plotly.react("cov_plot", c_traces, cov_layout)
    Plotly.react("scaled_plot", s_traces, scaled_layout)
}

\$(document).ready(function() {
    load_buttons(data.chromosomes)
    build_table()
})

</script>
</html>
"""


# variables from nextflow
bed = "$bedfile"
exclude = "$params.exclude"
ped = "$pedfile"
# known from indexcov
sample_col = "sample_id"
# known from indexcov
sex_col = "sex"
z_threshold = $params.zthreshold
distance_threshold = $params.distancethreshold
slop = $params.slop
gff = "$gff"
roc = "$rocfile"
# consistent between workflow and template
output = "covviz_report.html"

# start processing
logging.info("parsing bed file (%s)" % bed)
traces, samples = parse_bed(
    bed,
    exclude,
    ped,
    sample_col,
    sex_col,
    z_threshold,
    distance_threshold,
    slop,
)

logging.info("parsing gff file (%s)" % gff)
traces = parse_gff(gff, traces)
logging.info("parsing roc file (%s)" % roc)
traces = parse_roc(roc, traces, samples)
traces = parse_ped(ped, traces, sample_col)
with open(output, "w") as fh:
    html = TEMPLATE.replace("[DATA]", json.dumps(traces))
    print(html, file=fh)
logging.info("processing complete")
