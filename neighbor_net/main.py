# Run the neighbor net algorithm

import sys
import distances
import nnet_algorithm
import splits_io

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

    splits_io.print_splits_nexus(labels,splits,cycle)


if __name__ == '__main__':
    main()
