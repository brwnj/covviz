import csv
import logging
import os
import sys
import warnings
from collections import defaultdict
from itertools import groupby

import numpy as np
import pandas as pd

from .utils import gzopen

try:
    from itertools import ifilterfalse as filterfalse
except ImportError:
    from itertools import filterfalse

logger = logging.getLogger("covviz")


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
            logger.warning("%s is present in metadata, not in bed" % sample)
            valid = False
    if only_in_bed:
        for sample in only_in_bed:
            logger.warning("%s is present in bed, not in metadata" % sample)
            valid = False
    return valid


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
                        traces[sample]["x"].insert(
                            0, data["x"][index_values[0] - distance_idx]
                        )
                        traces[sample]["y"].insert(
                            0, data[sample][index_values[0] - distance_idx]
                        )
                    except IndexError as e:
                        # x_values[0] is the first data point
                        print(
                            "index error: %s\nindex_values[0]: %d, distance_idx: %d"
                            % (e, index_values[0], distance_idx)
                        )
                        break
                    extension_length -= (
                        data["x"][index_values[0] - distance_idx + 1]
                        - data["x"][index_values[0] - distance_idx]
                    )
                    distance_idx += 1

                traces[sample]["x"].extend(x_values)
                traces[sample]["y"].extend(y_values)

                # append slop
                extension_length = slop
                distance_idx = 1
                while extension_length > 0:
                    try:
                        traces[sample]["x"].append(
                            data["x"][index_values[-1] + distance_idx]
                        )
                        traces[sample]["y"].append(
                            data[sample][index_values[-1] + distance_idx]
                        )
                    except IndexError:
                        break
                    extension_length -= (
                        data["x"][index_values[-1] + distance_idx]
                        - data["x"][index_values[-1] + distance_idx - 1]
                    )
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
                logger.warning(
                    "sample column [%s] was not found in the header of %s"
                    % (sample_col, filename)
                )
                logger.warning(
                    "this will likely result in strange plotting behavior on sex chromosomes."
                )
                break
            if sex_col not in row:
                logger.warning(
                    "sex column [%s] was not found in the header of %s"
                    % (sex_col, filename)
                )
                logger.warning(
                    "this will likely result in strange plotting behavior on sex chromosomes."
                )
                break
            groups[row[sex_col]].append(row[sample_col])
    return groups


def normalize_depths(path, sex_chroms, median_window=7):
    filename, ext = os.path.splitext(path)
    if ext == ".gz":
        filename, ext = os.path.splitext(filename)

    output_bed = filename + ".norm.bed.gz"
    df = pd.read_csv(path, sep="\t", low_memory=False)
    # omit 0s from median calculation
    df[df.iloc[:, 3:] == 0] = np.nan
    # median values per sample
    global_sample_median = np.nanmedian(df.iloc[:, 3:], axis=0)
    # normalize within each sample
    df.iloc[:, 3:] = df.iloc[:, 3:] / global_sample_median

    extras = []
    for s in sex_chroms:
        if s.startswith("chr"):
            extras.append(s[3:])
        else:
            extras.append("chr" + s)
    sex_chroms = sex_chroms + extras

    # median per site
    autosome = ~np.asarray(df.iloc[:, 0].isin(sex_chroms))

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"All-NaN (slice|axis) encountered")
        site_median = np.nanmedian(df.iloc[:, 3:], axis=1)

    with np.errstate(invalid="ignore"):
        site_median[np.isnan(site_median) | (site_median < 0.03)] = 1
    # divide autosomes by median at each site so a given block is centered at
    # middle sample.
    for i in range(3, df.shape[1]):
        df.iloc[:, i] = np.where(autosome, df.iloc[:, i] / site_median, df.iloc[:, i])
        if median_window > 1:
            df.iloc[:, i] = (
                df.iloc[:, i].rolling(median_window).median()
            )  # pd.rolling_median(df.iloc[:, i], median_window)
        inan = np.asarray(np.isnan(df.iloc[:, i]))
        df.iloc[inan, i] = 0.0
    # df.to_csv(path_or_buf=output_bed, sep="\t", na_rep=0.0, index=False,
    #        compression='gzip',
    #        float_format="%.2f")
    return df
    # return output_bed


def identify_outliers(a, threshold=3.5):
    a = np.asarray(a, dtype=np.float32)
    med = np.median(a)
    mad = np.median(np.abs(a - med))
    # https://www.ibm.com/support/knowledgecenter/en/SSEP7J_11.1.0/com.ibm.swg.ba.cognos.ug_ca_dshb.doc/modified_z.html
    if mad == 0:
        meanAD = np.mean(np.abs(a - np.mean(a)))
        divisor = 1.253314 * meanAD
        modified_z_scores = (a - med) / divisor
    else:
        divisor = 1.4826 * mad
        modified_z_scores = (a - med) / divisor
    return np.where(np.abs(modified_z_scores) > threshold)


def add_roc_traces(path, traces, exclude):
    traces["roc"] = dict()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    n_bins = 150
    x_max = 2.5
    traces["roc"]["x_coords"] = [
        round(i, 2) for i in list(np.linspace(0, x_max, n_bins))
    ]

    for chrom, data in df.groupby(df.columns[0]):
        chrom = str(chrom)
        # apply exclusions
        if exclude.findall(chrom):
            continue

        chrom = chrom[3:] if chrom.startswith("chr") else chrom

        # pre-normalized data
        arr = np.asarray(data.iloc[:, 3:])
        traces["roc"][chrom] = dict()

        for i in range(0, arr.shape[1]):
            # get counts across our x-range of bins
            counts, _ = np.histogram(arr[:, i], bins=n_bins, range=(0, x_max))
            # decreasing order of the cumulative sum across the bins
            sums = counts[::-1].cumsum()[::-1]
            # normalize to y_max of 1
            sums = list(sums.astype(float) / max(1, sums[0]))
            traces["roc"][chrom][df.columns[i + 3]] = [round(i, 2) for i in sums]
    return traces


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
    min_samples=8,
    skip_norm=False,
    window=1,
):
    bed_traces = dict()
    # chromosomes, in order of appearance
    chroms = list()

    sex_chroms = [i.strip("chr") for i in sex_chroms.split(",")]

    groups = None
    if ped:
        groups = parse_sex_groups(ped, sample_col, sex_col)

    if skip_norm:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    else:
        df = normalize_depths(path, sex_chroms, median_window=window)

    if True:  # temporary to get sane diff.
        header = list(df.columns)
        samples = header[3:]
        if groups:
            valid = validate_samples(samples, groups)
            if not valid:
                logger.critical("sample ID mismatches exist between ped and bed")
                sys.exit(1)
        for chr, entries in df.groupby(header[0], as_index=False, sort=False):

            # apply exclusions
            if exclude.findall(chr):
                logger.debug("excluding chromosome: %s" % chr)
                continue

            data = defaultdict(list)
            bounds = dict(upper=[], lower=[])
            outliers = defaultdict(list)
            chrom = chr[3:] if chr.startswith("chr") else chr
            chroms.append(chrom)

            # capture plot area and outlier traces
            sample_groups = None
            if chrom in sex_chroms and groups:
                sample_groups = groups

            all_points = sample_groups is None

            if sample_groups is None:
                sample_groups = {"gid": samples}

            x_index = -1
            for not_x_index, row in entries.iterrows():
                x_index += 1

                x_value = int(row[header[1]])
                data["x"].append(x_value)

                for group_index, (gid, samples_of_group) in enumerate(
                    sample_groups.items()
                ):
                    # adds area traces where groups are present (chrs X and Y)
                    if len(bounds["upper"]) == group_index:
                        bounds["upper"].append([])
                        bounds["lower"].append([])
                    if all_points:
                        sample_values = np.minimum(
                            3, np.asarray(row[3:], dtype=np.float32)
                        )
                        for i, s in enumerate(samples):
                            data[s].append(float(sample_values[i]))
                    else:
                        sample_values = []
                        for sample in samples_of_group:
                            v = min(3.0, float(row[sample]))
                            data[sample].append(v)
                            sample_values.append(v)

                    # skip finding outliers for few samples
                    if len(samples) <= min_samples:
                        # save everything as an outlier
                        for sample in samples_of_group:
                            outliers[sample].append(
                                dict(index=x_index, x=x_value, y=data[sample][-1])
                            )
                        continue

                    # skip running test if everything is the same
                    if len(set(sample_values)) == 1:
                        bounds["upper"][group_index].append(float(sample_values[0]))
                        bounds["lower"][group_index].append(float(sample_values[0]))
                    else:
                        # indexes of passing values
                        passing = identify_outliers(sample_values, z_threshold)[0]
                        # remove those indexes from the list
                        # sometimes we get a list (x, y) for others we get a
                        # numpy array
                        if len(passing) > 0 and not isinstance(sample_values, list):
                            sample_values = sample_values.tolist()
                        for j in sorted(passing, reverse=True):
                            sample_values.pop(j)
                        # from remaining, grab upper and lower bounds
                        upper = float(max(sample_values))
                        lower = float(min(sample_values))
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
            json_output = dict(upper=[], lower=[], coords=data["x"], samples=[])

            # add the area traces
            for trace_index in range(len(bounds["upper"])):
                for bound in ["lower", "upper"]:
                    json_output[bound].append(
                        [round(i, 2) for i in bounds[bound][trace_index]]
                    )
            # add the sample traces for the outlier plots atop area traces
            len_traces = 0
            for sample, trace_data in traces.items():
                if not trace_data["x"]:
                    continue
                # y data may be gapped (string separated floats)
                y_data = list()
                for v in trace_data["y"]:
                    try:
                        y_data.append(round(v, 2))
                    except TypeError:
                        y_data.append(v)
                json_output["samples"].append(
                    {"name": sample, "x": trace_data["x"], "y": y_data}
                )
                len_traces += 1

            bed_traces[chrom] = json_output
            logger.info("plotting %d traces on chrom %s" % (len_traces, chr))

    bed_traces["chromosomes"] = chroms
    bed_traces["sample_list"] = samples
    # bed_traces["sex_chroms"] = sex_chroms

    # pass the bed or normed bed
    add_roc_traces(path, bed_traces, exclude)

    return bed_traces


def parse_bed_track(path, traces, exclude):
    """
    parse a bed file, placing lines per region
    """
    trace_name = os.path.basename(path).partition(".bed")[0]

    with gzopen(path) as fh:
        cleaned = filterfalse(lambda i: i[0] == "#", fh)

        for chrom, entries in groupby(
            cleaned, key=lambda i: i.partition("\t")[0].lstrip("chr")
        ):
            # apply exclusions
            if exclude.findall(chrom):
                continue
            if chrom not in traces:
                continue

            if not "annotations" in traces[chrom]:
                traces[chrom]["annotations"] = {"bed": []}
            if not "bed" in traces[chrom]["annotations"]:
                traces[chrom]["annotations"]["bed"] = list()

            regions = list()
            for line in entries:
                if line.startswith("#"):
                    continue
                toks = line.strip().split("\t")
                # not currently converting 0- and 1-based
                start = int(toks[1])
                end = int(toks[2])
                try:
                    name = toks[3]
                except IndexError:
                    name = ""
                regions.append([start, end, name])

            traces[chrom]["annotations"]["bed"].append([trace_name, regions])

    return traces
