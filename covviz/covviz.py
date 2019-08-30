"""
The aim of covviz is to highlight regions of significant (passing the
user's z-score threshold) and sustained (beyond user specified
distance) deviation from the majority of samples. Significance is
determined using z-scores for all samples at all points using median
absolute deviation, but in order to be highlighted, points must be
significant consecutively throughout a user specified distance.

If you are analyzing a low number of samples, deviation may be
irrelevant. In this case, we can set --min-samples to be greater
than our sample total to skip Z-threshold calculation and plot
coverages for all samples at all points.

Annotation tracks, --bed, --gff, and --vcf can be specified more than
once.
"""

import argparse
import csv
import json
import logging
import os
import re
from collections import defaultdict

from jinja2 import Environment, FileSystemLoader, select_autoescape
from lzstring import LZString

from .bed import parse_bed, parse_bed_track
from .gff import parse_gff
from .ped import parse_ped
from .utils import gzopen
from .vcf import parse_vcf

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("covviz")


class Formatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
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
    meta_group.add_argument(
        "--sex-vals",
        default="1,2",
        help="when using --ped, this defines male,female encoding",
    )

    annotations_group = p.add_argument_group("annotations")
    annotations_group.add_argument(
        "--bed",
        dest="bed_track",
        action="append",
        help=(
            ".bed reference file to place an optional region track; "
            "the name field is displayed on hover and overlapping "
            "regions will be merged."
        ),
    )
    annotations_group.add_argument(
        "--gff",
        action="append",
        help=(
            ".gff reference file to place an optional gene track; "
            "only rows of type 'gene' are used and annotated with "
            "gene ID where the attributes include 'Name=<symbol>;' "
            "annotation"
        ),
    )
    annotations_group.add_argument(
        "--gff-feature",
        default="gene",
        help="feature type within GFF upon which to search for attribute",
    )
    annotations_group.add_argument(
        "--gff-attr",
        default="Name=",
        help="the regex search string to grab the relevant portion of the attributes",
    )
    annotations_group.add_argument(
        "--vcf",
        action="append",
        help=(
            ".vcf file to place an optional variant track; info is "
            "displayed by default, but can be broken up using "
            "--vcf-info regex"
        ),
    )
    annotations_group.add_argument(
        "--vcf-info",
        default=None,
        help=(
            "the regex search string to grab the relevant portion of "
            "the INFO field, e.g. 'CLNDN=' for the clinical diagnosis "
            "field in ClinVar"
        ),
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

    annotation_track_colors = [
        "#666666",
        "#1b9e77",
        "#d95f02",
        "#7570b3",
        "#e7298a",
        "#66a61e",
        "#e6ab02",
        "#a6761d",
    ]

    y_offset = -0.15
    track_id = 0

    if args.gff:
        for gff in args.gff:
            logger.info("parsing gff file (%s)" % gff)
            traces = parse_gff(
                gff,
                traces,
                exclude,
                ftype=args.gff_feature,
                regex=args.gff_attr,
                y_offset=y_offset,
                track_color=annotation_track_colors[
                    track_id % len(annotation_track_colors)
                ],
            )
            y_offset += -0.15
            track_id += 1

    if args.bed_track:
        for bed in args.bed_track:
            logger.info("parsing bed file (%s)" % bed)
            traces = parse_bed_track(
                bed,
                traces,
                exclude,
                y_offset=y_offset,
                track_color=annotation_track_colors[
                    track_id % len(annotation_track_colors)
                ],
            )
            y_offset += -0.15
            track_id += 1

    if args.vcf:
        for vcf in args.vcf:
            logger.info("parsing vcf file (%s)" % vcf)
            traces = parse_vcf(
                vcf,
                traces,
                exclude,
                regex=args.vcf_info,
                y_offset=y_offset,
                track_color=annotation_track_colors[
                    track_id % len(annotation_track_colors)
                ],
            )
            y_offset += -0.15
            track_id += 1

    if args.ped:
        logger.info("parsing ped file (%s)" % args.ped)
        traces = parse_ped(
            args.ped, traces, args.sample_col, args.sex_chroms, args.sex_vals
        )

    with open(args.output, "w") as fh:
        logger.info("preparing output")
        compressed_json_data = compress_data(traces)
        html_template = env.get_template("covviz.html")
        print(html_template.render(data=compressed_json_data), file=fh)

    logger.info("processing complete")
