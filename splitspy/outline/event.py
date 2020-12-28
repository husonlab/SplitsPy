# event.py
"""An event used in the outline algorithm

See: Huson et al (2021)


LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""

from typing import Callable
from splitspy.splits.basic_split import Split

__author__ = "Daniel H. Huson"


class Event:
    def __init__(self, s: int, cycle: [int], splits: [Split], outbound: bool):
        self.__s = s
        self.__split = splits[s]
        self.__outbound = outbound
        self.__start_pos, self.__end_pos = splits[s].interval(cycle)

    def s(self) -> int:
        return self.__s

    def weight(self) -> float:
        return self.__split.weight

    def start_pos(self):
        return self.__start_pos

    def end_pos(self):
        return self.__end_pos

    def is_start(self):
        return self.__outbound

    def is_end(self):
        return not self.__outbound

    def __str__(self):
        return f'{self.__s} {self.__outbound}: {self.start_pos()} - {self.end_pos()} {self.weight(): .8f}'


def radix_sort(n_tax: int, outbound: [Event], inbound: [Event]) -> [Event]:
    outbound = __counting_sort(outbound, n_tax, lambda e: n_tax - e.end_pos())
    outbound = __counting_sort(outbound, n_tax, lambda e: e.start_pos())
    inbound = __counting_sort(inbound, n_tax, lambda e: n_tax - e.start_pos())
    inbound = __counting_sort(inbound, n_tax, lambda e: e.end_pos())

    return __merge(outbound, inbound)


def __counting_sort(events: [Event], max_key: int, key: Callable[[Event], int]) -> [Event]:
    if len(events) <= 1:
        return

    key2pos = [0] * (max_key + 1)

    for event in events:
        key2pos[key(event)] = key2pos[key(event)] + 1

    pos = 0
    for i in range(0, len(key2pos)):
        add = key2pos[i]
        key2pos[i] = pos
        pos += add

    other = [None] * len(events)

    for event in events:
        k = key(event)
        pos = key2pos[k]
        key2pos[k] = key2pos[k] + 1
        other[pos] = event

    return other


def __merge(outbound: [Event], inbound: [Event]) -> [Event]:
    ob = 0
    ib = 0

    events = []
    while ob < len(outbound) and ib < len(inbound):
        if outbound[ob].start_pos() < inbound[ib].end_pos() + 1:
            events.append(outbound[ob])
            ob += 1
        else:
            events.append(inbound[ib])
            ib += 1
    while ob < len(outbound):
        events.append(outbound[ob])
        ob += 1
    while ib < len(inbound):
        events.append(inbound[ib])
        ib += 1
    return events
