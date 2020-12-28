from typing import Tuple

__author__ = "David J. Bryant and Daniel H. Huson"

from splitpy.nnet import nnet_cycle, nnet_splits


def neighbor_net(labels: [str], mat: [float], cutoff=0.0001, constrained=True) -> Tuple[list, list]:
    cycle = nnet_cycle.compute(labels, mat)

    splits = nnet_splits.compute(len(labels), mat, cycle, cutoff, constrained)

    return cycle, splits
