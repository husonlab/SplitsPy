# nnet_alg.py
"""Runs the neighbor-net algorithm

See: Bryant and Moulton (2004)
See: Huson and Bryant (2006)


LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""
import time
from typing import Tuple

__author__ = "David J. Bryant and Daniel H. Huson"

from splitspy.nnet import nnet_cycle, nnet_splits, nnet_splits_lp


def neighbor_net(labels: [str], mat: [[float]], cutoff=0.0001, mode: str = "CLS", verbose: bool = True) -> Tuple[list, list]:
    if verbose:
        a = time.perf_counter()

    cycle = nnet_cycle.compute(labels, mat)

    if verbose:
        b = time.perf_counter()
        print(f"Computed circular order in {b-a:0.4f} seconds")

    if verbose:
        a = time.perf_counter()

    if mode == "LP":
        splits = nnet_splits_lp.compute(len(labels), mat, cycle, cutoff)
    else:
        constrained = (mode != "OLS")
        splits = nnet_splits.compute(len(labels), mat, cycle, cutoff, constrained)

    if verbose:
        b = time.perf_counter()
        print(f"Computed splits in {b-a:0.4f} seconds (mode={mode})")

    return cycle, splits
