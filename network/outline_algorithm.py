import math
from enum import Enum
from neighbor_net.csplit import CSplit
from network.graph import Graph
from network.event import Event
from network.event import radix_sort

__author__ = "Daniel H. Huson"


class Layout(Enum):
    UNROOTED = "unrooted"
    MIDPOINT = "midpoint"
    MIDPOINT_ALT = "midpoint_alt"
    OUTGROUP = "outgroup"
    OUTGROUP_ALT = "outgroup_alt"


def compute(labels: [str], cycle: [int], splits: [CSplit], layout=Layout.UNROOTED, out_group=None, used_wgts=True) -> Graph:
    n_tax = len(labels)

    splits = __add_trivial(n_tax, cycle, splits)

    split2angle = __compute_angles(n_tax, splits, 360.0 if layout == Layout.UNROOTED else 180.0)

    for s in range(0,len(splits)):
        print(splits[s]," -> ",split2angle[s])

    events = __setup_events(n_tax, splits)

    for e in events:
        print(e)

    xy = [0.0,0.0]

    graph = Graph()
    start = graph.new_node(location=xy,label=labels[0])

    current_splits = frozenset()

    splits2node = {current_splits: start}

    taxa_found = set()

    prev_event = None

    prev_node = start
    for event in events:
        if event.is_start():
            tmp = set(current_splits)
            tmp.add(event.s())
            current_splits=frozenset(tmp)
            xy = __translate(xy, split2angle[event.s()], event.weight() if used_wgts else 1)
        else:
            tmp = set(current_splits)
            tmp.remove(event.s())
            current_splits=frozenset(tmp)
            xy = __translate(xy, split2angle[event.s()] + 180, event.weight() if used_wgts else 1)

        need_new_node = (current_splits not in splits2node)

        if need_new_node:
            v = graph.new_node(location=xy)
            splits2node[current_splits] = v
        else:
            v = splits2node[current_splits]
            xy = v.location

        if not v.is_adjacent(prev_node):
            graph.new_edge(prev_node,v, None, event.weight() if used_wgts else 1,event.s())

            if not need_new_node:
                pass # just closed loop

        if prev_event is not None:
            if event.s() == prev_event.s():
                for t in splits[event.s()].part2():
                    v.label = labels[t-1]
                    taxa_found.add(t)

        prev_node = v
        prev_event = event

    return graph


def __add_trivial(n_tax: int, cycle: [int], splits: [CSplit]) -> [CSplit]:
    seen = set()

    for split in splits:
        if len(split.part1()) == 1:
            seen.add(i for i in split.part1())
        elif len(split.part2()) == 1:
            seen.add(i for i in split.part2())

    if len(seen) < n_tax:
        splits = splits.copy()
        for i in range(1,n_tax+1):
            if cycle[i] not in seen:
                splits.append(CSplit(n_tax,cycle,i,i,0))

    return splits


def __compute_angles(n_tax: int, splits: [CSplit], total_angle: float) -> [float]:
    angles = [0.0]
    for i in range(1,n_tax+1):
        angles.append((total_angle*(i-1.0)/n_tax)+270.0-0.5*total_angle)

    split_angle = []

    for s in splits:
        split_angle.append(__modulo360(0.5 * (angles[s.start_pos()] + angles[s.end_pos()])))

    return split_angle


def __modulo360(angle: float) -> float:
    while angle > 360:
        angle -= 360
    while angle < 0:
        angle += 360
    return angle


def __setup_events(n_tax: int, splits: [CSplit]):
    outbound = []
    inbound = []

    for s in range(0, len(splits)):
        outbound.append(Event(s, splits, True))
        inbound.append(Event(s, splits, False))

    return radix_sort(n_tax, outbound, inbound)

def __translate(xy: [float], angle: float, distance: float) -> [float]:
    dx = distance*math.cos( math.pi/180.0*angle)
    dy = distance*math.sin( math.pi/180.0*angle)

    if math.fabs(dx) < 0.000001:
        dx = 0.0
    if math.fabs(dy) < 0.000001:
        dy = 0.0
    return [xy[0]+dx,xy[1]+dy]









