# Run the neighbor net algorithm

import sys
import distances
import neighbor_net_algorithm

__author__ = "Daniel H. Huson"


def main():

    if len(sys.argv) != 3:
        print("Usage: neighbor_net_algorithm.py input outfile", file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    labels, matrix = distances.read(infile)

    distances.write(labels, matrix, outfile)

    splits, cycle = neighbor_net_algorithm.neighborNet(labels, matrix)

    print("Splits:")
    for split in splits:
        print(split)

    print("Cycle:",cycle)


if __name__ == '__main__':
    main()