from typing import Tuple

from split import Split
import neighbor_net_cycle
import neighbor_net_splits

__author__ = "Daniel H. Huson"


def neighbor_net(labels: [str], mat: [float], cutoff=0.0001, constrained=True) -> Tuple[list,list]:

    cycle = neighbor_net_cycle.compute(labels, mat)

    print("cycle:", cycle)

    splits = neighbor_net_splits.compute(len(labels), mat, cycle, cutoff,constrained)

    return splits, cycle
