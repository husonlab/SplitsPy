# outline_alg.py
"""Computes a phylogenetic outline from a set of circular splits

See: Bryant and Moulton (2004)
See: Huson et al (2021)


LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""
import math
from typing import Set, Tuple, List
from splitspy.splits import basic_split
from splitspy.splits.basic_split import Split, cyc_split
from splitspy.graph.graph import Graph
from splitspy.outlines.event import Event
from splitspy.outlines.event import radix_sort

__author__ = "Daniel H. Huson"


def compute(labels: [str], cycle: [int], splits: [Split], rooted: bool = False, alt: bool = False,
            out_grp: Set[int] = None, use_wts=True) -> Tuple[Graph, List[float]]:
    n_tax = len(labels)

    if rooted:
        if out_grp is None or len(out_grp) == 0:
            s, w1, w2 = __root_location_mid_point(alt, n_tax, cycle, splits, use_wts)
        else:
            s, w1, w2 = __root_location_out_group(out_grp, splits, use_wts)
        n_tax, labels, splits, cycle = __setup_rooted(alt, labels, splits, cycle, s, w1, w2)

    splits = __add_trivial(n_tax, cycle, splits)

    angles = __leaf_angles(n_tax, 360.0 if not rooted else 160.0)
    split2angle = __compute_angles(angles, cycle, splits)

    events = __setup_events(n_tax, cycle, splits)

    xy = [0.0, 0.0]

    graph = Graph()
    start = graph.new_node(pos=xy)

    current_splits = frozenset()

    splits2node = {current_splits: start}

    taxa_found = set()

    prev_event = None

    prev_node = start
    for event in events:
        if event.is_start():
            tmp = set(current_splits)
            tmp.add(event.s())
            current_splits = frozenset(tmp)
            xy = __translate(xy, split2angle[event.s()], event.weight() if use_wts else 1)
        else:
            tmp = set(current_splits)
            tmp.remove(event.s())
            current_splits = frozenset(tmp)
            xy = __translate(xy, split2angle[event.s()] + 180, event.weight() if use_wts else 1.0)

        need_new_node = (current_splits not in splits2node)

        if need_new_node:
            v = graph.new_node(pos=xy)
            splits2node[current_splits] = v
        else:
            v = splits2node[current_splits]
            xy = v.pos

        if not v.is_adjacent(prev_node):
            graph.new_edge(prev_node, v, None, event.weight() if use_wts else 1, event.s())

        if not need_new_node:
            pass  # just closed loop

        if prev_event is not None:
            if event.s() == prev_event.s():
                lab = []
                for t in splits[event.s()].part_not_in(cycle[1]):
                    lab.append(labels[t - 1])
                    taxa_found.add(t)
                prev_node.label = join(lab, ",")

        prev_node = v
        prev_event = event

    lab = []
    for t in range(1, n_tax+1):
        if t not in taxa_found:
            lab.append(labels[t-1])
    if len(lab) > 0:
        start.label = join(lab, ",")

    return graph, angles


def __add_trivial(n_tax: int, cycle: [int], splits: [Split]) -> [Split]:
    seen = set()

    for split in splits:
        if len(split.part1()) == 1:
            for t in split.part1():
                seen.add(t)
        elif len(split.part2()) == 1:
            for t in split.part2():
                seen.add(t)

    if len(seen) < n_tax:
        splits = splits.copy()
        for i in range(1, n_tax + 1):
            if cycle[i] not in seen:
                splits.append(cyc_split(cycle, i, i, 0))

    return splits


def __leaf_angles(n_tax: int, total_angle: float) -> [float]:
    angles = [0.0]
    for i in range(1, n_tax + 1):
        angles.append((total_angle * (i - 1.0) / n_tax) + 270.0 - 0.5 * total_angle)
    return angles


def __compute_angles(angles: [float], cycle: [int], splits: [Split]) -> [float]:

    split_angle = []

    for sp in splits:
        a, b = sp.interval(cycle)
        split_angle.append(__modulo360(0.5 * (angles[a] + angles[b])))

    return split_angle


def __modulo360(angle: float) -> float:
    while angle >= 360.0:
        angle -= 360.0
    while angle < 0.0:
        angle += 360.0
    return angle


def __setup_events(n_tax: int, cycle: [int], splits: [Split]):
    outbound = []
    inbound = []

    for s in range(0, len(splits)):
        outbound.append(Event(s, cycle, splits, True))
        inbound.append(Event(s, cycle, splits, False))

    return radix_sort(n_tax, outbound, inbound)


def __translate(xy: [float], angle: float, distance: float) -> [float]:
    dx = distance * math.cos(math.pi / 180.0 * angle)
    dy = distance * math.sin(math.pi / 180.0 * angle)

    if math.fabs(dx) < 0.000001:
        dx = 0.0
    if math.fabs(dy) < 0.000001:
        dy = 0.0
    return [xy[0] + dx, xy[1] + dy]


def __root_location_mid_point(alt: bool, n_tax: int, cycle: [int], splits: [Split], use_wts: bool)\
        -> Tuple[int, float, float]:
    dist = basic_split.split_dist(n_tax, splits)

    max_dist = 0.0
    furthest = [0, 0]

    for a in range(1,n_tax+1):
        for b in range(a+1,n_tax+1):
            if dist[a][b] > max_dist:
                max_dist = dist[a][b]
                furthest = [min(a, b), max(a, b)]

    split2idx = {}

    sep = []

    for s in range(0, len(splits)):
        sp = splits[s]
        if sp.separates(furthest[0], furthest[1]):
            sep.append(sp)
            split2idx[sp] = s

    interval = __interval(furthest[0], furthest[1], cycle, alt)

    triples = []
    pos = 0
    for sp in sep:
        p = sp.part_in(furthest[0])
        triples.append([len(interval.intersection(p)), len(p), split2idx[sp]])
        pos += 1

    triples.sort()

    total = 0.0
    for trp in triples:
        sp = splits[trp[2]]
        wgt = sp.weight if use_wts else 1.0
        delta = total + wgt - 0.5 * max_dist
        if delta > 0:
            return split2idx[sp], delta, wgt - delta
        total += wgt
    return 1, 0.0, splits[0].weight if use_wts else 1.0


def __interval(a: int, b: int, cycle: [int], alt: bool) -> Set[int]:
    interval = set()

    if len(cycle) > 0:
        if alt:
            entered = False
            i = len(cycle) - 1
            while True:
                if cycle[i] == a:
                    interval.add(a)
                    entered = True
                if entered and cycle[i] == b:
                    break
                if i == 1:
                    i = len(cycle) - 1
                else:
                    i -= 1
        else:
            entered = False
            i = 1
            while True:
                if cycle[i] == a:
                    interval.add(a)
                    entered = True
                if entered and cycle[i] == b:
                    break
                if i >= len(cycle) - 1:
                    i = 1
                else:
                    i += 1
    return interval


def __root_location_out_group(out_grp: Set[int],  splits: [Split],use_wts: bool) -> Tuple[int, float, float]:

    if len(out_grp) > 0:
        out_grp_splits = set()
        out_grp_taxon = min(out_grp)

        for p in range(0, len(splits)):
            pa = splits[p].part_in(out_grp_taxon)
            if out_grp <= pa:
                ok = True
                to_delete = set()
                for q in out_grp_splits:
                    qa = splits[q].part_in(out_grp_taxon)
                    if qa <= pa:
                        ok = False
                        break
                    elif pa <= qa:
                        to_delete.add(q)
                if len(to_delete) > 0:
                    out_grp_splits.difference_update(to_delete)
                if ok:
                    out_grp_splits.add(p)

        if len(out_grp_splits) > 0:
            s = min(out_grp_splits)
            return s, 0.9*splits[s].weight, 0.1*splits[s].weight

    return 1, 0.0, splits[0].weight if use_wts else 1.0


def __setup_rooted(alt: bool, labels0: [str], splits0: [Split], cycle0: [int], mid: int, w1: float, w2: float) -> Tuple:
    labels = labels0.copy()
    labels.append("Root")
    n_tax = len(labels)
    root_id = n_tax

    cycle = [0] * (len(cycle0) + 1)

    first = 0
    if not alt:
        part = splits0[mid].part_not_in(1)
        t = 1
        for v in cycle0:
            if v > 0:
                if first == 0 and v in part:
                    first = v
                    cycle[t] = root_id
                    t = t + 1
                cycle[t] = v
                t += 1
    else:
        part = splits0[mid].part_not_in(1)
        seen = 0
        t = 1
        for v in cycle0:
            if v > 0:
                cycle[t] = v
                t += 1
                if v in part:
                    seen += 1
                    if seen == len(part):
                        first = v
                        cycle[t] = root_id
                        t += 1

    cycle = rotate(cycle, root_id)

    total_wgt = 0.0

    mid1 = splits0[mid].deepcopy()
    mid1.part_in(1).add(root_id)
    mid1.weight = w1

    mid2 = splits0[mid].deepcopy()
    mid2.part_not_in(1).add(root_id)
    mid2.weight = w2

    splits = []
    for s in range(0, len(splits0)):
        if s == mid:
            total_wgt += mid1.weight
            splits.append(mid1)
        else:
            sp = splits0[s].deepcopy()
            p = mid1.part_not_in(root_id)
            if set(sp.part1()) <= p:
                sp.part2().add(root_id)
            elif set(sp.part2()) <= p:
                sp.part1().add(root_id)
            elif len(sp.part_in(first)) > 1:
                sp.part_in(first).add(root_id)
            else:
                sp.part_not_in(first).add(root_id)
            splits.append(sp)
            total_wgt += sp.weight

    total_wgt += mid2.weight
    splits.append(mid2)

    splits.append(cyc_split(cycle, 2, n_tax, total_wgt / len(splits) if total_wgt > 0 else 1.0))

    return n_tax, labels, splits, cycle


def rotate(cycle: [int], first: int) -> [int]:
    result = [0]
    for i in range(1, len(cycle)):
        if cycle[i] == first:
            while len(result) < len(cycle):
                result.append(cycle[i])
                i += 1
                if i == len(cycle):
                    i = 1
            break
    return result


def join(lst: [str], sep: str):
    s = ""
    for a in lst:
        if len(s) > 0:
            s += sep
        s += a
    return s


