# Run the neighbor net algorithm

import sys
import phylo_outline.utils.distances as distances
import phylo_outline.nnet.nnet_algorithm as nnet_algorithm
from phylo_outline.utils import splits_io
import phylo_outline.outline.outline_algorithm
from phylo_outline.utils.graph_utils import networkx

__author__ = "Daniel H. Huson"


def main():
    if len(sys.argv) != 3:
        print("Usage: nnet_algorithm.py input outfile", file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    labels, matrix = distances.read(infile)

    # distances.write(labels, matrix, outfile)

    cycle, splits = nnet_algorithm.neighbor_net(labels, matrix)

    splits_io.print_splits_nexus(labels, splits, cycle)

    graph = phylo_outline.outline.outline_algorithm.compute(labels, cycle, splits)

    graph.write_tgf()

    nx_graph = networkx(graph)

if __name__ == '__main__':
    main()
