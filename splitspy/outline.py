# outline.py
"""Program that computes a phylogenetic outline for a distance matrix, using neighbor-net

Given a distance matrix on a set of taxa, this program runs the neighbor-net algorithm
and then computes a phylogenetic outline. The outline is drawn in a graphics window

See: Bryant and Moulton (2004)
See: Huson et al (2021)


LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""

import splitspy.nnet.distances as distances
import splitspy.nnet.nnet_algo as nnet_algorithm
from splitspy.graph import draw
from splitspy.splits import splits_io
import splitspy.outline.outline_algo
from optparse import OptionParser, OptionGroup
from splitspy.splits.basic_split import split_dist

__author__ = "Daniel H. Huson"


def main():
    parser = OptionParser("%prog [options] infile",
                          description="Run neighbor-net and compute a phylogenetic outline",
                          epilog="Please cite: Huson et al (2021) and Bryant and Moulton (2004)."
                                 " Drawing uses John Zelle's graphics package.")

    parser.add_option("-o", "--output", default="", action="store", dest="outfile", help="Output file",
                      metavar="FILE")

    parser.add_option("-n", "--no_draw", default=False, action="store_true", dest="no_draw",
                      help="Do not draw the network")

    outline_opts = OptionGroup(parser, "Outline Options")
    outline_opts.add_option("-r", "--rooted", default=False, action="store_true", dest="rooted", help="Rooted network")

    outline_opts.add_option("-a", "--alt", default=False, action="store_true", dest="alt",
                            help="Alternative layout rooted network")

    outline_opts.add_option("-g", "--ogroup", default="", action="store", type="string", dest="out_grp",
                            help="Out-group taxa for rooted network (format: tax1,tax2,...)", metavar="GRP")

    parser.add_option_group(outline_opts)

    win_opts = OptionGroup(parser, "Window Options")
    win_opts.add_option("--width", default=1000, action="store", dest="win_width", help="Window width",
                        metavar="WIDTH", type="int")
    win_opts.add_option("--height", default=800, action="store", dest="win_height", help="Window height",
                        metavar="HEIGHT", type="int")
    win_opts.add_option("--m_left", default=100, action="store", dest="ml", help="Left margin", type="int",
                        metavar="MARGIN")
    win_opts.add_option("--m_right", default=100, action="store", dest="mr", help="Right margin", type="int",
                        metavar="MARGIN")
    win_opts.add_option("--m_top", default=100, action="store", dest="mt", help="Top margin", type="int",
                        metavar="MARGIN")
    win_opts.add_option("--m_bot", default=100, action="store", dest="mb", help="Bottom margin", type="int",
                        metavar="MARGIN")
    win_opts.add_option("--font_size", default=12, action="store", dest="f_size", help="Font size", type="int")

    parser.add_option_group(win_opts)

    (options, args) = parser.parse_args()

    if len(args) == 1:
        infile = args[0]
    else:
        infile = "-"

    if len(args) >= 2:
        raise IOError("Too many arguments", args)

    if options.outfile == "":
        outfile = None
    else:
        outfile = options.outfile

    labels, matrix = distances.read(infile)

    # distances.write(labels, matrix, outfile)

    cycle, splits = nnet_algorithm.neighbor_net(labels, matrix)

    fit = distances.ls_fit(matrix, split_dist(len(labels), splits))

    if outfile is not None:
        splits_io.print_splits_nexus(labels, splits, cycle, fit, filename=outfile)

    if options.out_grp != "":
        names = set(options.out_grp.split(","))
        unknown = names.difference(set(labels))
        if len(unknown) > 0:
            raise IOError("Unknown taxa in out-group:", unknown)
        out_grp = set()
        for t in range(1, len(labels)+1):
            if labels[t-1] in names:
                out_grp.add(t)
    else:
        out_grp = None

    graph, angles = splitspy.outline.outline_algo.compute(labels, cycle, splits, rooted=options.rooted,
                                                          out_grp=out_grp, alt=options.alt)

    # graph.write_tgf()

    if not options.no_draw:
        title = infile if infile != "-" else "Phylogenetic outline"
        draw.draw(graph, angles, fit, title, options.win_width, options.win_height, options.ml,
                  options.mr, options.mt, options.mb, options.f_size)


if __name__ == '__main__':
    main()
