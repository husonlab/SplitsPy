from typing import Tuple

import nnet_cycle
import nnet_splits

__author__ = "Daniel H. Huson"


def neighbor_net(labels: [str], mat: [float], cutoff=0.0001, constrained=True) -> Tuple[list,list]:

    cycle = nnet_cycle.compute(labels, mat)

    splits = nnet_splits.compute(len(labels), mat, cycle, cutoff, constrained)

    return cycle, splits
