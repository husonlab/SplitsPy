# nnet_splits_lp.py
"""Computes splits weights using LP

See: Bryant and Huson, 2021


LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""
from typing import List

from scipy.optimize import linprog
import numpy as np
from splitspy.splits.basic_split import *


__author__ = "Daniel H. Huson"


def compute(n_tax: int, mat: np.array, cycle: [int], cutoff=0.00001) -> [Split]:
    """ compute splits and their weights using Linear Program
        Parameters
        ----------
            n_tax: int
                number of taxa
            mat: np.array
                distance matrix, 0-based
            cycle: [int]
                circular ordering, 1-based
            cutoff: float
                minimum split weight
        Returns
        -------
            [Split]
                splits with weights, taxa are 1-based
          """

    if n_tax == 1:
        return []
    elif n_tax == 2:
        return [cyc_split([1, 2], 2, 2, mat[0][1])] if mat[0][1] >= cutoff else []

    all = __all_splits(cycle)

    total = 0.0
    split_counts = np.zeros(len(all),dtype=np.int32)
    for i in range(0, n_tax):
        for j in range(i + 1, n_tax):
            total += mat[i][j]
            for s in range(0, len(all)):
                split = all[s]
                if split.separates(i+1, j+1):
                    split_counts[s] += 1

    c = []

    for s in range(0, len(all)):
       c.append(-split_counts[s])

    A = []
    b = []

    for i in range(0, n_tax):
        for j in range(i + 1, n_tax):
            row = []
            for s in range(0, len(all)):
                split = all[s]
                if split.separates(i + 1, j + 1):
                    row.append(1)
                else:
                    row.append(0)
            A.append(row)
            b.append(mat[i][j])

    res = linprog(c, A_ub=A, b_ub=b)

    result = []
    for s in range(0, len(all)):
        # print("w"+str(s+1)+":",res.x[s])
        if res.x[s] > cutoff:
            src = all[s]
            split = Split(src.part1(),src.part2(),res.x[s])
            result.append(split)

    print(f"Delta: {total+res.fun:0.4f}")

    return result


def __all_splits(cycle: List[int]) -> List[Split]:
    splits = []

    for p in range(2, len(cycle)):
        for q in range(p, len(cycle)):
            splits.append(cyc_split(cycle,p,q,1))

    return splits


