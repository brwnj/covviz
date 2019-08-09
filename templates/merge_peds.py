#!/usr/bin/env python
import csv
import gzip
import logging

custom_ped_file = "$custom_ped"
standard_ped_file = "$ped"
output_file = "merged.ped"
sample_col = "$params.samplecol"
omit_from_indexcov = ["#family_id", "paternal_id", "maternal_id", "phenotype", "p.out", "PC4", "PC5"]
sep = "\\t"

def gzopen(f):
    if f.endswith(".gz"):
        return gzip.open(f, "rt")
    else:
        return open(f)


indexcov_sample_col = "sample_id"
custom_header = []
custom_data = dict()
table_data = []

with gzopen(custom_ped_file) as fh:
    custom_header = fh.readline().strip().split(sep)
    if "sample_id" in custom_header:
        omit_from_indexcov.append("sample_id")
    fh.seek(0)
    reader = csv.DictReader(fh, delimiter=sep)
    for row in reader:
        try:
            custom_data[row[sample_col]] = row
        except KeyError:
            logging.critical(" you may need to change `--samplecol` to reflect your sample ID column")
            raise

with gzopen(standard_ped_file) as fh:
    header = fh.readline().strip().split(sep)
    header = [i for i in header if i not in omit_from_indexcov]
    merged_header = custom_header + header
    table_data.append(merged_header)
    fh.seek(0)
    reader = csv.DictReader(fh, delimiter=sep)
    for row in reader:
        merged_row = []

        if row[indexcov_sample_col] not in custom_data:
            logging.warning("sample %s was not present in %s" % (row[indexcov_sample_col], custom_ped_file))
            # output will need to be padded on the left
            for col in merged_header:
                try:
                    merged_row.append(row[col])
                except KeyError:
                    merged_row.append("")
        else:
            # custom data
            sample_data = custom_data[row[indexcov_sample_col]]
            # indexcov data
            for col in header:
                sample_data[col] = row[col]
            for col in merged_header:
                merged_row.append(sample_data[col])

        table_data.append(merged_row)

with open(output_file, "w") as fh:
    for line in table_data:
        print(*line, sep=sep, file=fh)
