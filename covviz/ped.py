import csv
from collections import defaultdict

from .utils import gzopen


def check_ped_header(header, cols):
    """
    header (list) - current ped header columns
    cols (list) - typical columns we'll see from indexcov
    """
    header = set(header)
    cols = set(cols)
    is_indexcov = True
    if len(cols - header) > 0:
        is_indexcov = False
    return is_indexcov


def parse_ped(path, traces, sample_col, sex_chroms, sex_vals="1,2"):
    table_data = list()
    ped_data = dict(
        inferred=defaultdict(list), bins=defaultdict(list), pca=defaultdict(list)
    )
    sex_chroms = [i.strip() for i in sex_chroms.split(",")]

    male = None
    female = None
    for i, v in enumerate(sex_vals.split(",")):
        if i == 0:
            male = v
        elif i == 1:
            female = v
        else:
            break

    # indexcov .ped format
    required_vals = [
        sample_col,
        "bins.in",
        "bins.out",
        "bins.lo",
    ]  # CNX, CNY or CNchrX, CNchrY
    float_vals = ["slope", "p.out", "PC1", "PC2", "PC3", "PC4", "PC5"]
    int_vals = ["bins.out", "bins.lo", "bins.hi", "bins.in"]

    via_indexcov = False

    with gzopen(path) as fh:
        header = fh.readline().strip().split("\t")
        # datatables needs these escaped
        datatable_cols = [dict(title=i, data=i.replace(".", "\\.")) for i in header]
        table_data.append(datatable_cols)
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter="\t")

        via_indexcov = check_ped_header(header, required_vals)

        for row in reader:

            # coerce vals
            if via_indexcov:
                for k, v in row.items():
                    if k in int_vals:
                        row[k] = int(v)
                    elif k in float_vals:
                        row[k] = float(v)
                    elif k.startswith("CN"):
                        row[k] = float(v)

            # coerced vals saved in data.ped
            table_data.append(row)

            if via_indexcov:
                # inferred sex
                ped_data["inferred"]["x"].append(row["CN%s" % sex_chroms[0]])
                try:
                    ped_data["inferred"]["y"].append(row["CN%s" % sex_chroms[1]])
                except IndexError:
                    ped_data["inferred"]["y"].append(0)

                # male
                if row["sex"] == male:
                    ped_data["inferred"]["color"].append("rgba(12,44,132,0.5)")
                    ped_data["inferred"]["hover"].append(
                        "Sample: %s<br>Inferred X CN: 1" % (row[sample_col],)
                    )
                # female
                elif row["sex"] == female:
                    ped_data["inferred"]["color"].append("rgba(227,26,28,0.5)")
                    ped_data["inferred"]["hover"].append(
                        "Sample: %s<br>Inferred X CN: 2" % (row[sample_col],)
                    )
                else:
                    ped_data["inferred"]["color"].append("rgba(105,105,105,0.5)")
                    ped_data["inferred"]["hover"].append(
                        "Sample: %s<br>Unknown" % (row[sample_col],)
                    )

                # bin plot
                total = row["bins.in"] + row["bins.out"]
                ped_data["bins"]["samples"].append(row[sample_col])
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
    if via_indexcov:
        traces["depth"] = ped_data
    return traces
