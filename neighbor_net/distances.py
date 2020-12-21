import sys
from typing import Tuple

import numpy as np

__author__ = 'Daniel Huson'


def read(filename="-") -> Tuple[list,np.array]:
    if filename == "-":
        ins = sys.stdin
    else:
        ins = open(filename)

    n = 0

    labels = []
    rows = []

    for line in ins:
        if n == 0:
            n = int(line)
            if n <= 0:
                raise IOError("Number of taxa must be positive, got: " + line)
        else:
            tokens = line.split()
            if len(tokens) != n + 1:
                raise IOError("Wrong number of tokens in line, got: " + line)
            labels.append(tokens[0])
            row = []
            for i in range(1, n + 1):
                row.append(float(tokens[i]))
            rows.append(row)

    if ins != sys.stdin:
        ins.close()

    return labels, np.array(rows)


def write(labels: [str], matrix: [float], filename="-") -> None:
    if filename == "-":
        outs = sys.stdout
    else:
        outs = open(filename, mode="w")
    n = len(labels)
    print(n, file=outs)

    for i in range(0, n):
        print(labels[i], end="", file=outs)
        row = matrix[i]
        for value in row:
            print(" ", value, end="", file=outs)
        print(end="\n", file=outs)

    if outs != sys.stdout:
        outs.close()
