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
    sex_chroms="X,Y",
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
    sex_chroms = [i.strip("chr") for i in sex_chroms.split(",")]

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
    bed_traces["sex_chroms"] = sex_chroms
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


def parse_ped(path, traces, sample_col, sex_chroms):
    table_data = list()
    ped_data = dict(inferred=defaultdict(list), bins=defaultdict(list), pca=defaultdict(list))
    sex_chroms = [i.strip() for i in sex_chroms.split(",")]
    float_vals = ["slope", "p.out", "PC1", "PC2", "PC3", "PC4", "PC5"]
    int_vals = ["sex", "bins.out", "bins.lo", "bins.hi", "bins.in"]

    with gzopen(path) as fh:
        header = fh.readline().strip().split("\t")
        datatable_cols = [dict(title=i, data=i.replace(".", "\\.")) for i in header]
        table_data.append(datatable_cols)
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:

            # coerce vals
            for k, v in row.items():
                if k in int_vals:
                    row[k] = int(v)
                elif k in float_vals:
                    row[k] = float(v)
                elif k.startswith("CN"):
                    row[k] = float(v)

            # coerced vals saved in data.ped
            table_data.append(row)
            # inferred sex
            ped_data["inferred"]["x"].append(row["CN%s" % sex_chroms[0]])
            try:
                ped_data["inferred"]["y"].append(row["CN%s" % sex_chroms[1]])
            except IndexError:
                ped_data["inferred"]["y"].append(0)
            # male
            if row["sex"] == 1:
                ped_data["inferred"]["color"].append('rgba(12,44,132,0.5)')
                ped_data["inferred"]["hover"].append('Sample: %s<br>Inferred X CN: 1' % (row["sample_id"],))
            else:
                ped_data["inferred"]["color"].append('rgba(227,26,28,0.5)')
                ped_data["inferred"]["hover"].append('Sample: %s<br>Inferred X CN: 2' % (row["sample_id"],))
            # bin plot
            total = row["bins.in"] + row["bins.out"]
            ped_data["bins"]["samples"].append(row["sample_id"])
            ped_data["bins"]["x"].append(row["bins.lo"] / total)
            ped_data["bins"]["y"].append(row["bins.out"] / total)
            # PCAs
            try:
                ped_data["pca"]["pca_1"].append(row["PC1"])
                ped_data["pca"]["pca_2"].append(row["PC2"])
                ped_data["pca"]["pca_3"].append(row["PC3"])
            except KeyError:
                pass

    traces["ped"] = table_data
    traces["sample_column"] = sample_col
    traces["depth"] = ped_data
    return traces


TEMPLATE = """<!DOCTYPE html>
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
    .nav-pills .nav-link.active, .nav-pills .show > .nav-link {background-color: #6c757d !important;}
    table.dataTable thead th.sorting:after,
    table.dataTable thead th.sorting_asc:after,
    table.dataTable thead th.sorting_desc:after,
    table.dataTable thead th.sorting:before,
    table.dataTable thead th.sorting_asc:before,
    table.dataTable thead th.sorting_desc:before {
        font-family: FontAwesome !important;
    }
    </style>
</head>
<body>
    <nav class="navbar navbar-dark navbar-expand-md p-0" style="background-color: #3C444C;">
        <a class="navbar-brand m-0 p-1 text-light"><img class="pr-4" src='data:image/png;base64,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' alt="" height="58">
        </a>
        <div class="collapse navbar-collapse">
            <ul class="nav nav-pills" role="tablist">
                <li class="nav-item">
                    <a class="nav-link text-white active" data-toggle="tab" role="tab" href="#chrom_coverage" aria-controls="chrom_coverage" aria-selected="true">Chromosome Coverage</a>
                </li>
                <li class="nav-item">
                    <a class="nav-link text-white" data-toggle="tab" role="tab" href="#global_qc" aria-controls="global_qc" aria-selected="false">Global QC</a>
                </li>
            </ul>
        </div>
        <span class="navbar-brand mb-0 h1">covviz report</span>
    </nav>
    <div class="container-fluid w-90 pt-3">
        <div class="tab-content border-bottom">
            <div class="tab-pane fade show active" id="chrom_coverage" role="tabpanel" aria-labelledby="scaled_tab">
                <div class="row">
                    <div class="col-12">
                        <div id="chrom_selector" class="btn-group btn-group-toggle d-flex justify-content-center flex-wrap" data-toggle="buttons">
                            <div class="input-group-prepend">
                                <span class="input-group-text">Chromosome:</span>
                            </div>
                        </div>
                    </div>
                </div>
                <ul class="nav nav-tabs pt-3" role="tablist">
                    <li class="nav-item">
                        <a class="nav-link active text-secondary" id="scaled_tab" data-toggle="tab" href="#scaled" role="tab" aria-controls="scaled" aria-selected="true">Scaled chromosome coverage</a>
                    </li>
                    <li class="nav-item">
                        <a class="nav-link text-secondary" id="cov_tab" data-toggle="tab" href="#cov" role="tab" aria-controls="cov" aria-selected="false">Proportions covered</a>
                    </li>
                </ul>
                <div class="tab-content border-bottom">
                    <div class="tab-pane fade show active mb-3" id="scaled" role="tabpanel" aria-labelledby="scaled_tab">
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
                    <div class="tab-pane fade mb-3" id="cov" role="tabpanel" aria-labelledby="cov_tab">
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
            </div>
            <div class="tab-pane fade" id="global_qc" role="tabpanel" aria-labelledby="global_qc">
                <div id="global_qc_plots" class="container-fluid w-90">
                    <div class="row border-bottom">
                        <div class="col-6">
                            <h4>Inferred sex</h4>
                        </div>
                        <div class="col-6">
                            <h4>Problematic low and non-uniform coverage bins</h4>
                        </div>
                    </div>
                    <div class="row">
                        <div id="inferred_sex" class="col-6"></div>
                        <div id="bin_counts" class="col-6"></div>
                    </div>
                </div>
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
var scatter_point_color = 'rgba(31,120,180,0.5)'
var scatter_point_color_light = 'rgba(255,255,255,0.1)'

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
    cov_plot.removeAllListeners("plotly_doubleclick")
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
    scaled_plot.removeAllListeners("plotly_doubleclick")
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
                    reset_line_plots()
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

const get_sex_data = (x, y, text, hover, colors) => {
    return {
        x: x,
        y: y,
        mode: 'markers',
        type: 'scatter',
        text: text,
        hovertext: hover,
        marker: { size: 12, color: colors, line: { color: 'rgb(40,40,40)', width: 1, } },
    }
}

const get_bin_data = (x, y, text, color='rgba(31,120,180,0.5)') => {
    return {
        x: x,
        y: y,
        mode: 'markers',
        type: 'scatter',
        text: text,
        hoverinfo: 'text',
        marker: { size: 12, color: color, line: {color: 'rgb(40,40,40)', width: 1,} },
    }
}

const get_pca_data = (x, y, text, color='rgba(31,120,180,0.5)') => {
    return {
        x: x,
        y: y,
        mode: 'markers',
        type: 'scatter',
        text: text,
        hoverinfo: 'text',
        marker: { size: 12, color: color, line: { color: 'rgb(40,40,40)', width: 1, } },
    }
}

const build_inferred_sex_plot = () => {
    var sex_p = document.getElementById("inferred_sex")
    Plotly.react(sex_p,
        [get_sex_data(data.depth.inferred.x, data.depth.inferred.y, data.depth.bins.samples, data.depth.inferred.hover, data.depth.inferred.color)],
        {
            margin: {t: 15,},
            height: 400,
            xaxis: {title: "X Copy Number"},
            yaxis: {title: "Y Copy Number"},
            hovermode: "closest",
            legend: { orientation: "h", x: 0, y: 1.125, },
            dragmode: "lasso",
        }
    )
    sex_p.removeAllListeners("plotly_selected")
    sex_p.removeAllListeners("plotly_deselect")
    sex_p.on("plotly_selected", handle_scatter_selection)
    sex_p.on("plotly_deselect", handle_scatter_deselect)
}

const build_global_qc = () => {

    // inferred sex
    build_inferred_sex_plot()

    // bin counts
    var bin_p = document.getElementById("bin_counts")
    Plotly.react(bin_p, [get_bin_data(data.depth.bins.x, data.depth.bins.y, data.depth.bins.samples)], {
        margin: {t: 15,},
        height: 400,
        xaxis: {title: "Proportion of bins with depth < 0.15"},
        yaxis: {title: "Proportion of bins with depth<br>outside of (0.85, 1.15)"},
        hovermode: "closest",
        legend: { orientation: "h", x: 0, y: 1.125, },
        dragmode: "lasso",
    })
    bin_p.on("plotly_selected", handle_scatter_selection)
    bin_p.on("plotly_deselect", handle_scatter_deselect)

    // PCAs
    if ("pca_1" in data.depth.pca) {
        // insert HTML into plot area
        let pa = document.getElementById("global_qc_plots")
        pa.insertAdjacentHTML('beforeend', '<div class="row border-bottom"><div class="col-6"><h4>PCA 1 vs 2</h4></div><div class="col-6"><h4>PCA 1 vs 3</h4></div></div><div class="row"><div id="pca_1" class="col-6"></div><div id="pca_2" class="col-6"></div></div>')
        let pca1_p = document.getElementById("pca_1")
          pca2_p = document.getElementById("pca_2")
        Plotly.react(pca1_p, [get_pca_data(data.depth.pca.pca_1, data.depth.pca.pca_2, data.depth.bins.samples)], {
            margin: {t: 15,},
            height: 400,
            xaxis: {title: "PC1"},
            yaxis: {title: "PC2"},
            hovermode: "closest",
            legend: { orientation: "h", x: 0, y: 1.125, },
            dragmode: "lasso",
        })
        pca1_p.on("plotly_selected", handle_scatter_selection)
        pca1_p.on("plotly_deselect", handle_scatter_deselect)

        Plotly.react(pca2_p, [get_pca_data(data.depth.pca.pca_1, data.depth.pca.pca_3, data.depth.bins.samples)], {
            margin: {t: 15,},
            height: 400,
            xaxis: {title: "PC1"},
            yaxis: {title: "PC3"},
            hovermode: "closest",
            legend: { orientation: "h", x: 0, y: 1.125, },
            dragmode: "lasso",
        })
        pca2_p.on("plotly_selected", handle_scatter_selection)
        pca2_p.on("plotly_deselect", handle_scatter_deselect)
    }
}

\$('#chrom_selector').on("change", () => {
    // \$('#tab_names a[href="#scaled"]').tab('show')
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

const search_datatable = (sample_id) => {
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

const reset_datatable = () => {
    let tables = \$('.dataTable').DataTable()
    let table = tables.table('#ped_table')
    table.\$('tr.selected').removeClass('selected')
    table.search('').columns().search('').draw()
}

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
    search_datatable(sample_id)
            }
        }
    }
}

const handle_scatter_deselect = (event) => {
    build_inferred_sex_plot()
    Plotly.restyle('bin_counts', 'marker.color', scatter_point_color, [0])
    if ("pca_1" in data.depth.pca) {
        Plotly.restyle('pca_1', 'marker.color', scatter_point_color, [0])
        Plotly.restyle('pca_2', 'marker.color', scatter_point_color, [0])
    }
}

const handle_scatter_selection = (event) => {
    var sample_ids = event.points.map(point => point.text)
    var colors = []
    for (var i = 0; i < data.depth.bins.samples.length; i++) {
        colors.push(scatter_point_color_light)
    }
    event.points.forEach((pt) => {
        colors[pt.pointNumber] = scatter_point_color
    })

    Plotly.restyle('inferred_sex', 'marker.color', [colors])
    Plotly.restyle('bin_counts', 'marker.color', [colors])
    if ("pca_1" in data.depth.pca) {
        Plotly.restyle('pca_1', 'marker.color', [colors])
        Plotly.restyle('pca_2', 'marker.color', [colors])
    }
}

const reset_line_plots = () => {
    scaled_layout.xaxis.autorange = true
    // de-select scaled plot
    Plotly.react("scaled_plot", scaled_traces, scaled_layout)
    cov_layout.xaxis.range = [0, 1.5]
    cov_layout.yaxis.range = [0, 1.]
    // de-select cov plot
    Plotly.react("cov_plot", cov_traces, cov_layout)
}

const handle_plot_doubleclick = () => {
    reset_line_plots()
    // de-select in table
    if ("ped" in data) {
        reset_datatable()
    }
}

const highlight_plot_traces = (sample_id) => {
    let s_traces = []
    let c_traces = []
    let k_traces = []
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
    build_global_qc()
})

</script>
</html>
"""

# variables from nextflow
bed = "$bedfile"
exclude = "$params.exclude".replace("~", "").replace(",", "|")
ped = "$pedfile"
# known from indexcov
sample_col = "sample_id"
# known from indexcov
sex_col = "sex"
sex_chroms = "$params.sexchroms"
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
    sex_chroms,
    z_threshold,
    distance_threshold,
    slop,
)
logging.info("parsing gff file (%s)" % gff)
traces = parse_gff(gff, traces)
logging.info("parsing roc file (%s)" % roc)
traces = parse_roc(roc, traces, samples)
traces = parse_ped(ped, traces, sample_col, sex_chroms)
with open(output, "w") as fh:
    html = TEMPLATE.replace("[DATA]", json.dumps(traces))
    print(html, file=fh)
logging.info("processing complete")

