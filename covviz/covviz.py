import argparse
import csv
import json
import logging
import os
from collections import defaultdict

from jinja2 import Environment, FileSystemLoader, select_autoescape
from lzstring import LZString

from .bed import parse_bed
from .gff import parse_gff
from .ped import parse_ped
from .utils import gzopen

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("covviz")


# TODO calculate this based on the bed
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
                marker={"color": "rgba(108,117,125,0.2)"},
            )
            traces["roc"][chr].append(trace)
    return traces


def compress_data(data):
    json_str = json.dumps(data).encode("utf-8", "ignore").decode("utf-8")
    json_str = json_str.replace("NaN", "null")
    return LZString().compressToBase64(json_str)


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # this will later be the only requirement
    p.add_argument("bed", help="indexcov .bed output")

    p.add_argument(
        "--exclude",
        default="^GL|^hs|^chrEBV$|M$|MT$|^NC|_random$|Un_|^HLA\-|_alt$|hap\d+$",
        help="chromosome regex to exclude from analysis",
    )
    p.add_argument("--ped", help="indexcov .ped output")
    p.add_argument(
        "--sex-chroms", default="X,Y", help="sex chromosomes, e.g. chrX,chrY"
    )
    p.add_argument(
        "--z-threshold",
        default=3.5,
        type=float,
        help="the point at which we determine a sample is an outlier from the group at any given point",
    )
    p.add_argument(
        "--distance-threshold",
        default=150000,
        type=int,
        help="when calculating significance, the Z-score has to be above the threshold for consecutive points up to the total distance set by distance threshold",
    )
    p.add_argument(
        "--slop",
        default=500000,
        type=int,
        help="slop is the distance to add to traces when plotting -- without slop, it's not always clear what happens to the points immediately flanking the area of significant deviation",
    )
    p.add_argument("--gff", help=".gff reference file")
    p.add_argument("--roc", help="indexcov .roc output")
    p.add_argument("--output", default="covviz_report.html", help="output file path")
    p.add_argument(
        "--min-samples",
        default=6,
        type=int,
        help="show all traces when analyzing this few samples; ignores z-threshold, distance-threshold, and slop",
    )
    return p.parse_args()


def cli():
    args = parse_args()


    env = Environment(
        loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')), autoescape=select_autoescape(["html"])
    )

    exclude = args.exclude.replace("~", "").replace(",", "|")
    sample_col = "sample_id"
    sex_col = "sex"

    logger.info("parsing bed file (%s)" % args.bed)
    traces, samples = parse_bed(
        args.bed,
        args.exclude,
        args.ped,
        sample_col,
        sex_col,
        args.sex_chroms,
        args.z_threshold,
        args.distance_threshold,
        args.slop,
        args.min_samples,
    )
    logger.info("parsing gff file (%s)" % args.gff)
    traces = parse_gff(args.gff, traces)
    logger.info("parsing roc file (%s)" % args.roc)
    traces = parse_roc(args.roc, traces, samples)
    traces = parse_ped(args.ped, traces, sample_col, args.sex_chroms)
    with open(args.output, "w") as fh:
        compressed_json_data = compress_data(traces)
        html_template = env.get_template("covviz.html")
        print(html_template.render(data=compressed_json_data), file=fh)
    logger.info("processing complete")
