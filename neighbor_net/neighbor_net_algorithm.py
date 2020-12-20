from split import Split
import neighbor_net_cycle

__author__ ="Daniel H. Huson"


def neighborNet (labels, matrix):

    cycle = neighbor_net_cycle.compute_neighbor_net_cycle(labels,matrix)

    splits = []
    n = len(labels)
    for i in range(1,n+1):
        splits.append(Split(range(1,i),range(i,n+1),1))

    return splits, cycle