import gzip


def gzopen(f):
    if f.endswith(".gz"):
        return gzip.open(f, "rt")
    else:
        return open(f)
