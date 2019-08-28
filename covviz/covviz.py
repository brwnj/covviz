import argparse
import csv
import json
import logging
import os
import re
from collections import defaultdict

from jinja2 import Environment, FileSystemLoader, select_autoescape
from lzstring import LZString

from .bed import parse_bed
from .gff import parse_gff
from .ped import parse_ped
from .utils import gzopen

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("covviz")


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # this will later be the only requirement
    p.add_argument(
        "bed",
        help=(
            "bed3+ file format with a header defining sample columns "
            "after chrom, start, and end."
        ),
    )
    p.add_argument(
        "-e",
        "--exclude",
        default="^HLA,^hs,:,^GL,M,EBV,^NC,^phix,decoy,random$,Un,hap,_alt$",
        help="chromosome regex to exclude from analysis",
    )
    p.add_argument(
        "-x",
        "--sex-chroms",
        default="X,Y",
        help="sex chromosomes as they are defined in your bed, e.g. chrX,chrY",
    )

    p.add_argument(
        "-z",
        "--z-threshold",
        default=3.5,
        type=float,
        help=(
            "the point at which we determine a sample is an outlier "
            "from the group at any given point"
        ),
    )
    p.add_argument(
        "-d",
        "--distance-threshold",
        default=150000,
        type=int,
        help=(
            "when calculating significance, the Z-score has to be "
            "above the threshold for consecutive points up to the "
            "total distance set by distance threshold"
        ),
    )
    p.add_argument(
        "-s",
        "--slop",
        default=500000,
        type=int,
        help=(
            "slop is the distance to add to traces when plotting -- "
            "without slop, it's not always clear what happens to the "
            "points immediately flanking the area of significant "
            "deviation"
        ),
    )
    p.add_argument(
        "--gff",
        help=(
            ".gff reference file to place an optional gene track; "
            "only rows of type 'gene' are used and annotated with "
            "gene ID where the attributes include 'Name=<symbol>;' "
            "annotation"
        ),
    )

    p.add_argument(
        "-o", "--output", default="covviz_report.html", help="output file path"
    )
    p.add_argument(
        "--skip-norm",
        action="store_true",
        help=(
            "skip normalization by global sample median if the depths "
            "in your .bed are already normalized"
        ),
    )
    p.add_argument(
        "--min-samples",
        default=8,
        type=int,
        help=(
            "show all traces when analyzing this few samples; ignores "
            "z-threshold, distance-threshold, and slop"
        ),
    )

    meta_group = p.add_argument_group("sample metadata")
    meta_group.add_argument(
        "-p", "--ped", help="ped file defining samples, sex, and other metadata"
    )
    meta_group.add_argument(
        "--sample-col",
        default="sample_id",
        help="when using --ped, this defines the sample ID column",
    )
    meta_group.add_argument(
        "--sex-col", default="sex", help="when using --ped, this defines the sex column"
    )

    return p.parse_args()


def compress_data(data):
    json_str = json.dumps(data).encode("utf-8", "ignore").decode("utf-8")
    json_str = json_str.replace("NaN", "null")
    return LZString().compressToBase64(json_str)


def cli():
    args = parse_args()

    env = Environment(
        loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), "templates")),
        autoescape=select_autoescape(["html"]),
    )

    logger.info("parsing bed file (%s)" % args.bed)

    exclude = re.compile(args.exclude.replace("~", "").replace(",", "|"))

    traces, samples = parse_bed(
        args.bed,
        exclude,
        args.ped,
        args.sample_col,
        args.sex_col,
        args.sex_chroms,
        args.z_threshold,
        args.distance_threshold,
        args.slop,
        args.min_samples,
        args.skip_norm,
    )

    if args.gff:
        logger.info("parsing gff file (%s)" % args.gff)
        traces = parse_gff(args.gff, traces, exclude)

    if args.ped:
        logger.info("parsing ped file (%s)" % args.ped)
        traces = parse_ped(args.ped, traces, args.sample_col, args.sex_chroms)

    with open(args.output, "w") as fh:
        logger.info("preparing output")
        compressed_json_data = compress_data(traces)
        html_template = env.get_template("covviz.html")
        print(html_template.render(data=compressed_json_data), file=fh)

    logger.info("processing complete")
