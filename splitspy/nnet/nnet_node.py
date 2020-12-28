# nnet_node.py
"""A node used by the neighbor-net algorithm

See: Bryant and Moulton (2004)
See: Huson and Bryant (2006)


LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""
__author__ = "David J. Bryant and Daniel H. Huson"


class NetNode:
    def __init__(self, _id):
        self.id = _id
        self.nbr = None
        self.ch1 = None
        self.ch2 = None
        self.next = None
        self.prev = None
        self.Rx = 0.0
        self.Sx = 0.0

    def create_string(self) -> str:
        string = "[id=" + str(self.id) + ", nbr="
        string += "null" if self.nbr is None else str(self.nbr.id)
        string += ", ch1="
        string += "null" if self.ch1 is None else str(self.ch1.id)
        string += ", ch2="
        string += "null" if self.ch2 is None else str(self.ch2.id)
        string += ", prev="
        string += "null" if self.prev is None else str(self.prev.id)
        string += ", next="
        string += "null" if self.next is None else str(self.next.id)
        string += ", Rx=" + str(self.Rx)
        string += ", Sx=" + str(self.Sx)
        string += "]"
        return string

    def __repr__(self):
        return self.create_string()
