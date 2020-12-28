# splits_io.py
"""Basic output of splits

LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""

import sys
import splitspy
from splitspy.splits import basic_split

__author__ = "Daniel H. Huson"


def print_splits_nexus(labels: [str], splits: [basic_split.Split], cycle: [int], fit=-1.0, filename="-") -> None:
    if filename == "-":
        outs = sys.stdout
    else:
        outs = open(filename, mode="w")

    print("#nexus", file=outs)

    print("BEGIN taxa;", file=outs)
    print("DIMENSIONS nTax=", len(labels), ";", sep="", file=outs)
    print("TAXLABELS", file=outs)
    for label in labels:
        print("'", label, "'", sep="", file=outs)
    print(";", file=outs)
    print("END;", file=outs)

    print("BEGIN SPLITS;", file=outs)
    print("DIMENSIONS nTax=", len(labels), " nSplits=", len(splits), ";", sep="", file=outs)
    print("FORMAT labels=no weights=yes confidences=no;", file=outs)
    print("PROPERTIES", end=" ", file=outs)
    if fit != -1:
        print("fit=", fit, end=" ", file=outs)
    print("compatible" if splitspy.splits.basic_split.compatible(splits) else "cyclic", ",", sep="", file=outs)
    print("CYCLE", end="", file=outs)
    for i in range(1, len(cycle)):
        print(" ", cycle[i], sep="", end="", file=outs)
    print(";", file=outs)
    print("MATRIX", file=outs)
    for sp in splits:
        print(f'{sp.get_weight():.8f}', end="\t", file=outs)
        first = True
        for t in sp.part1():
            if first:
                first = False
            else:
                print(" ", end="", file=outs)
            print(t, end="", file=outs)
        print(",", file=outs)
    print(";", file=outs)
    print("END;", file=outs)

    if outs != sys.stdout:
        outs.close()


def print_splits_fasta(labels: [str], splits: [basic_split.Split], filename="-") -> None:
    if filename == "-":
        outs = sys.stdout
    else:
        outs = open(filename, mode="w")

    for i in range(0, len(labels)):
        print(">", labels[i], sep="", file=outs)
        for split in splits:
            print("1" if (i + 1) in split.part1() else "0", end="", file=outs)
        print(file=outs)

    if outs != sys.stdout:
        outs.close()
