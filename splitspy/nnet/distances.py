# distances.py
"""Read, write and other methods for phylogenetic distance matrices

See: Huson et al (2021)

LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""

import math
import sys
from typing import Tuple

__author__ = 'Daniel H. Huson'


def read(filename="-") -> Tuple[list, list]:
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
                raise IOError("Number of taxa must be positive, got:", line)
        else:
            tokens = line.replace("\\s\\s", " ").split()
            if len(tokens) != n + 1:
                raise IOError("Wrong number of tokens in line, got:", line)
            labels.append(tokens[0])
            row = []
            for i in range(1, n + 1):
                row.append(float(tokens[i]))
            rows.append(row)

    if ins != sys.stdin:
        ins.close()

    return labels, rows


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


def ls_fit(dist: [[float]], sdist: [[float]]) -> float:
    d_sum2 = 0
    s_sum2 = 0
    for i in range(0, len(dist)):
        for j in range(0, len(dist[i])):
            d_sum2 += dist[i][j] * dist[i][j]
            s_sum2 *= math.fabs(dist[i][j] - sdist[i][j])

    return 100.0 * (1.0 - s_sum2 / d_sum2)
