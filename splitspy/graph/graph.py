# graph.py
"""A simple graph implementation using linked lists

See: Huson et al (2021)


LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""

import sys

__author__ = 'Daniel H. Huson'


class Node:
    def __init__(self, _id: int, label: str = None, pos: [float] = None, info=None):
        self.__id = _id
        self.label = label
        self.pos = pos
        self.info = info

        self.__prev = None
        self.__next = None

        self.__first_src_edge = None
        self.__first_tar_edge = None

        self.__in_deg = 0
        self.__out_deg = 0

    def __str__(self):
        return f'{self.__id} {self.label}' if self.pos is None else f'{self.__id} {self.label} {self.pos[0]} {self.pos[1]}'

    def deg(self) -> int:
        return self.__in_deg + self.__out_deg

    def in_deg(self) -> int:
        return self.__in_deg

    def out_deg(self) -> int:
        return self.__out_deg

    def in_edges(self) -> []:
        e = self.__first_tar_edge
        while e is not None:
            yield e
            e = e._Edge__next_tar_edge

    def out_edges(self):
        e = self.__first_src_edge
        while e is not None:
            yield e
            e = e._Edge__next_src_edge

    def adj_edges(self) -> []:
        e = self.__first_tar_edge
        while e is not None:
            yield e
            e = e._Edge__next_tar_edge
        e = self.__first_src_edge
        while e is not None:
            yield e
            e = e._Edge__next_src_edge

    def is_adjacent(self, other) -> bool:
        for e in self.adj_edges():
            if e.src() == other or e.tar() == other:
                return True
        return False

    def children(self):
        e = self.__first_tar_edge
        while e is not None:
            yield e.tar()
            e = e._Edge__next_tar_edge

    def parents(self) -> []:
        e = self.__first_tar_edge
        while e is not None:
            yield e.src()
            e = e._Edge__next_tar_edge

    def id(self):
        return self.__id

    def next(self):
        return self.__next

    def prev(self):
        return self.__prev

    def __hash__(self) -> int:
        return hash(self.__id)

    def __eq__(self, other):
        return self.__id == other.__id


class Edge:
    def __init__(self, _id: int, src: Node, tar: Node, label: str = None, weight: float = 1, info=None):
        self.__id = _id
        self.__src = src
        self.__tar = tar
        self.label = label
        self.weight = weight
        self.info = info

        self.__prev = None
        self.__next = None

        self.__next_src_edge = None
        self.__prev_src_edge = None
        self.__next_tar_edge = None
        self.__prev_tar_edge = None

    def __str__(self):
        return f'{self.__src.id()} {self.__tar.id()} {self.weight}'

    def src(self) -> Node:
        return self.__src

    def tar(self) -> Node:
        return self.__tar

    def opp(self, v: Node) -> Node:
        if self.__src is v:
            return self.__tar
        elif self.__tar is v:
            return self.__src
        else:
            return v

    def id(self):
        return self.__id

    def next(self):
        return self.__next

    def prev(self):
        return self.__prev

    def print_info(self):
        print(self, end="")
        if self.__prev is not None:
            print(" prev =", self.__prev.id(), end="")
        if self.__next is not None:
            print(" next =", self.__next.id(), end="")
        print()

    def __hash__(self) -> int:
        return hash(self.__id)

    def __eq__(self, other):
        return self.__id == other.__id


class Graph(object):
    def __init__(self):
        self.__first_node = None
        self.__last_node = None
        self.__first_edge = None
        self.__last_edge = None

        self.__top_node_id = 0
        self.__top_edge_id = 0

        self.__n_nodes = 0
        self.__n_edges = 0

    def n_nodes(self) -> int:
        return self.__n_nodes

    def n_edges(self) -> int:
        return self.__n_edges

    def clear(self) -> None:
        self.__first_node = None
        self.__last_node = None
        self.__first_edge = None
        self.__last_edge = None

        self.__top_node_id = 0
        self.__top_edge_id = 0

        self.__n_nodes = 0
        self.__n_edges = 0

    def new_node(self, label: str = None, pos: [float] = None, info=None) -> Node:
        self.__top_node_id += 1
        v = Node(self.__top_node_id, label, pos, info)

        if self.__first_node is None:
            self.__first_node = v
            self.__last_node = v
        else:
            if self.__last_node is not None:
                v._Node__prev = self.__last_node
                self.__last_node._Node__next = v
            self.__last_node = v
        self.__n_nodes += 1

        return v

    def new_edge(self, src: Node, tar: Node, label: str = None, weight: float = -1, info=None) -> Edge:
        self.__top_edge_id += 1
        edge = Edge(self.__top_edge_id, src, tar, label, weight, info)

        if self.__first_edge is None:
            self.__first_edge = edge
            self.__last_edge = edge
        else:
            edge._Edge__prev = self.__last_edge
            self.__last_edge._Edge__next = edge
            self.__last_edge = edge

        p = src._Node__first_src_edge
        if p is None:
            src._Node__first_src_edge = edge
        else:
            while p._Edge__next_src_edge is not None:
                p = p._Edge__next_src_edge
            p._Edge__next_src_edge = edge
            edge._Edge__prev_src_edge = p

        q = tar._Node__first_tar_edge
        if q is None:
            tar._Node__first_tar_edge = edge
        else:
            while q._Edge__next_tar_edge is not None:
                q = q._Edge__next_tar_edge
            q._Edge__next_tar_edge = edge
            edge._Edge__prev_tar_edge = q

        self.__n_edges += 1

        return edge

    def delete_edge(self, edge: Edge) -> None:
        src = edge.src()
        if src._Node__first_src_edge is edge:
            src._Node__first_src_edge = edge._Edge__next_src_edge
        else:
            edge._Edge__prev_src_edge._Edge__next_src_edge = edge._Edge__next_src_edge

        if edge._Edge__next_src_edge is not None:
            edge._Edge__next_src_edge._Edge__prev_src_edge = edge._Edge__prev_src_edge

        tar = edge.tar()
        if tar._Node__first_tar_edge is edge:
            tar._Node__first_tar_edge = edge._Edge__next_tar_edge
        else:
            edge._Edge__prev_tar_edge._Edge__next = edge._Edge__next_tar_edge

        if edge._Edge__next_tar_edge is not None:
            edge._Edge__next_tar_edge._Edge__prev = edge._Edge__prev_tar_edge

        if self.__first_edge is edge:
            self.__first_edge = edge._Edge__next
        if self.__last_edge is edge:
            self.__last_edge = edge._Edge__prev

        if edge._Edge__prev is not None:
            h = edge._Edge__next
            k = edge._Edge__prev
            k._Edge__next = h

        if edge._Edge__next is not None:
            edge._Edge__next._Edge__prev = edge._Edge__prev

        self.__n_edges -= 1

    def delete_node(self, node: Node) -> None:
        for e in node.adj_edges():
            self.delete_edge(e)

        if self.__first_node is node:
            self.__first_node = node._Node__next
        if self.__last_node is node:
            self.__last_node = node._Node__prev

        if node._Node__prev is not None:
            node._Node__prev._Node__next = node._Node__next

        if node._Node__next is not None:
            node._Node__next._Node__prev = node._Node__prev

        self.__n_nodes -= 1

    def first_node(self):
        return self.__first_node

    def last_node(self):
        return self.__last_node

    def first_edge(self):
        return self.__first_edge

    def last_edge(self):
        return self.__last_edge

    def nodes(self):
        v = self.__first_node
        while v is not None:
            yield v
            v = v._Node__next

    def edges(self):
        e = self.__first_edge
        while e is not None:
            yield e
            e = e._Edge__next

    def bbox(self):
        x_min = 1000000.0
        x_max = -1000000.0
        y_min = 1000000.0
        y_max = -1000000.0

        for v in self.nodes():
            x_min = min(x_min, v.pos[0])
            x_max = max(x_max, v.pos[0])
            y_min = min(y_min, v.pos[1])
            y_max = max(y_max, v.pos[1])

        return x_min, x_max, y_min, y_max

    def write_tgf(self, outfile="-") -> None:
        if outfile == "-":
            outs = sys.stdout
        else:
            outs = open(outfile, mode="w")

        for v in self.nodes():
            print(v.id(), end="", file=outs)
            if v.label is not None:
                print(" ", v.label, end="", file=outs)
            if v.pos is not None:
                print(" [", "{:.6f}".format(v.pos[0]), ",", "{:.6f}".format(v.pos[1]), "]",
                      sep="", end="", file=outs)
            if v.info is not None:
                print(" {", v.info, "}", sep="", end="", file=outs)
            print()

        for e in self.edges():
            print(e.src().id(), e.tar().id(), end="", file=outs)
            if e.weight != -1:
                print(" [", "{:.6f}".format(e.weight), "]", sep="", end="", file=outs)
            if e.info is not None:
                print(" {", e.info, "}", sep="", end="", file=outs)
            print()

        if outs != sys.stdout:
            outs.close()


if __name__ == "__main__":
    g = Graph()

    nodes = [None]
    weight = 1
    for i in range(1, 5):
        nodes.append(g.new_node("x" + str(i), [i, i]))
        for j in range(1, i):
            g.new_edge(nodes[j], nodes[i], weight=weight)
            weight += 1

    print("nodes", g.n_nodes(), "edges", g.n_edges())
    g.write_tgf()

    for v in g.nodes():
        print(v, "in:")
        for e in v.in_edges():
            print("\t", e)

    for v in g.nodes():
        print(v, "out:")
        for e in v.out_edges():
            print("\t", e)

    g.delete_edge(list(g.edges())[0])

    print("nodes", g.n_nodes(), "edges", g.n_edges())
    g.write_tgf()

    g.delete_edge(list(g.edges())[2])

    print("nodes", g.n_nodes(), "edges", g.n_edges())
    g.write_tgf()

    g.delete_node(list(g.nodes())[3])

    print("nodes", g.n_nodes(), "edges", g.n_edges())
    g.write_tgf()
