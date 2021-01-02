# outline.py
"""Program that computes a phylogenetic outline for a distance matrix, using neighbor-net

Given a distance matrix on a set of taxa, this program runs the neighbor-net algorithm
and then computes a phylogenetic outline. The outline is drawn in a graphics window

See: Bryant and Moulton (2004)
See: Huson et al (2021)


LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""
from typing import Set

import splitspy.nnet.distances as distances
import splitspy.nnet.nnet_algo as nnet_algorithm
from splitspy.graph import draw
from splitspy.splits import splits_io
import splitspy.outlines.outline_algo
from optparse import OptionParser, OptionGroup
from splitspy.splits.basic_split import split_dist

__author__ = "Daniel H. Huson"


def main():
    """ Run neighbor net and compute a phylogenetic outline

    Usage:
    -----
    python splitspy.outline.py [options] infile

    Options:
    -------
        -h, --help          show this help message and exit
        -o FILE, --output=FILE
                            output image file
        -n FILE, --nexus=FILE
                            output splits file (Nexus format for SplitsTree5)
        -t FILE, --tgf=FILE output graph file (in trivial graph format)

        Outline Options:
        -r, --rooted        rooted network
        -a, --alt           alternative layout for rooted network
        -g GRP, --out_grp=GRP
                            out-group taxa for rooted network (format: tax1,tax2,...)

        Window Options:
        --width=WIDTH       window width
        --height=HEIGHT     window height
        --m_left=MARGIN     left margin
        --m_right=MARGIN    right margin
        --m_top=MARGIN      top margin
        --m_bot=MARGIN      bottom margin
        --font_size=F_SIZE  font size

    Input format:
    ------------
    Example:

        6
        A.andrenof  0 0.090103 0.103397 0.096012 0.004431 0.075332
        A.mellifer  0.090103 0 0.093058 0.090103 0.093058 0.100443
        A.dorsata   0.103397 0.093058 0 0.116691 0.106352 0.103397
        A.cerana    0.096012 0.090103 0.116691 0 0.098966 0.09896
        A.florea    0.004431 0.093058 0.106352 0.098966 0 0.078287
        A.koschev   0.075332 0.100443 0.103397 0.098966 0.078287 0

    Contributions:
    -------------
    The neighbor-net algorithm is due to David J. Bryant and Vincent Moulton (2004).
    It was originally implemented in Matlab by David Bryant. David Bryant and Daniel Huson ported the code to Java
    (Huson and Bryant, 2006) and to Python (Huson et al, 2021).
    Phylogenetic outlines and the outline algorithm are due to David Bryant and Daniel Huson, and were implemented
    by Daniel Huson in Java and Python (Huson et al, 2021).
"""
    parser = OptionParser("%prog [options] infile",
                          description="Run neighbor-net and compute a phylogenetic outline",
                          epilog="Please cite: Huson et al (2021) and Bryant and Moulton (2004).")

    parser.add_option("-o", "--output", default="", action="store", dest="outfile", help="output image file",
                      metavar="FILE")

    parser.add_option("-n", "--nexus", default="", action="store", dest="nexus_file",
                      help="output splits file (Nexus format for SplitsTree5)",  metavar="FILE")

    parser.add_option("-t", "--tgf", default="", action="store", dest="graph_file",
                      help="output graph file (in trivial graph format)",  metavar="FILE")

    outline_opts = OptionGroup(parser, "Outline Options")
    outline_opts.add_option("-r", "--rooted", default=False, action="store_true", dest="rooted", help="rooted network")

    outline_opts.add_option("-a", "--alt", default=False, action="store_true", dest="alt",
                            help="alternative layout for rooted network")

    outline_opts.add_option("-g", "--out_grp", default="", action="store", type="string", dest="out_grp_labels",
                            help="out-group taxa for rooted network (format: tax1,tax2,...)", metavar="GRP")

    parser.add_option_group(outline_opts)

    win_opts = OptionGroup(parser, "Window Options")
    win_opts.add_option("--width", default=1000, action="store", dest="win_width", help="window width",
                        metavar="WIDTH", type="int")
    win_opts.add_option("--height", default=800, action="store", dest="win_height", help="window height",
                        metavar="HEIGHT", type="int")
    win_opts.add_option("--m_left", default=100, action="store", dest="m_left", help="left margin", type="int",
                        metavar="MARGIN")
    win_opts.add_option("--m_right", default=100, action="store", dest="m_right", help="right margin", type="int",
                        metavar="MARGIN")
    win_opts.add_option("--m_top", default=100, action="store", dest="m_top", help="top margin", type="int",
                        metavar="MARGIN")
    win_opts.add_option("--m_bot", default=100, action="store", dest="m_bot", help="bottom margin", type="int",
                        metavar="MARGIN")
    win_opts.add_option("--font_size", default=12, action="store", dest="font_size", help="font size", type="int",
                        metavar="SIZE")

    parser.add_option_group(win_opts)

    (options, args) = parser.parse_args()

    if len(args) == 1:
        infile = args[0]
    elif len(args) == 0:
        raise IOError("Must specify input file (use - for stdin)")
    else:
        raise IOError("Too many arguments", args)

    labels, matrix = distances.read(infile)

    out_grp = set()
    if options.out_grp_labels is not None and options.out_grp_labels != "":
        out_grp_labels = set(options.out_grp_labels.split(","))
        unknown = out_grp_labels.difference(set(labels))
        if len(unknown) > 0:
            raise IOError("Unknown taxa in out-group:", unknown)
        for t in range(1, len(labels) + 1):
            if labels[t - 1] in out_grp_labels:
                out_grp.add(t)

    run(labels, matrix, outfile=options.outfile, nexus_file=options.nexus_file, graph_file=options.graph_file,
        rooted=options.rooted, alt=options.alt, out_grp=out_grp,
        win_width=options.win_width,win_height=options.win_height, m_left=options.m_left, m_right=options.m_right,
        m_top=options.m_top, m_bot=options.m_bot,font_size=options.font_size)


def run(labels: [str], matrix: [[float]], outfile: str = "", nexus_file: str = "", graph_file: str = "",
        rooted: bool = False, alt: bool = False, out_grp: Set[int] = None, win_width: int = 1000, win_height: int = 800,
        m_left: int = 100, m_right: int = 100, m_top: int = 100, m_bot: int = 100, font_size: int = 12) -> None:


    # distances.write(labels, matrix, outfile)

    cycle, splits = nnet_algorithm.neighbor_net(labels, matrix)

    fit = distances.ls_fit(matrix, split_dist(len(labels), splits))

    if nexus_file != "":
        splits_io.print_splits_nexus(labels, splits, cycle, fit, filename=nexus_file)

    graph, angles = splitspy.outlines.outline_algo.compute(labels, cycle, splits, rooted=rooted, out_grp=out_grp, alt=alt)

    if graph_file != "":
        graph.write_tgf(outfile=graph_file)

    draw.draw(outfile, graph, angles, fit, win_width, win_height,m_left, m_right, m_top, m_bot, font_size)


if __name__ == '__main__':
    main()
