#!/usr/bin/env python
# coding=utf-8

import argparse


def main(html, py):
    html_str = ""
    py_str = ""
    with open(html) as fh:
        for line in fh:
            line = line.replace("$", "\$")
            html_str += line
    with open(py) as fh:
        for line in fh:
            if line.strip() == "TEMPLATE = [TEMPLATE]":
                line = 'TEMPLATE = """%s"""' % (html_str,)
            py_str += line
    with open("parse_indexcov.py", "w") as fh:
        print(py_str, file=fh)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("html", help="covviz.html component")
    p.add_argument("py", help="covviz.py component")
    args = p.parse_args()
    main(args.html, args.py)
