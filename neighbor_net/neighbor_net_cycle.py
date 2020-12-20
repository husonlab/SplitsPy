import numpy as np
from collections import deque

__author__ = 'Daniel Huson'


def compute_neighbor_net_cycle(labels, matrix):
    n = len(labels)

    if n <= 3:
        return list(range(0, n + 1))

    nodes_head = setup_nodes(n)

    mat = setup_matrix(labels, matrix) # matrix is 0-based, mat is 1-based

    # print("Mat:", mat)

    joins = join_nodes(n, mat, nodes_head)

    cycle = expand_nodes(joins, nodes_head)

    return cycle


def setup_nodes(n):
    nodes_head = NetNode(0)

    for i in range(n, 0, -1):
        node = NetNode(i)
        node.next = nodes_head.next
        nodes_head.next = node

    node = nodes_head
    while node.next is not None:
        node.next.prev = node
        node = node.next

    return nodes_head


def setup_matrix(labels, matrix):
    n = len(labels)
    max_number_of_nodes = max(3, 3 * n - 5)

    values = []

    for i in range(0, max_number_of_nodes+1):
        values.append(0.0)

    for i in range(0, max_number_of_nodes):
        values.append(0.0)
        for j in range(0, max_number_of_nodes):
            if i < n and j < n:
                values.append(matrix[i][j])
            else:
                values.append(0)

    tmp = labels.copy()
    for i in range(n, max_number_of_nodes):
        tmp.append(str(i))

    return np.array(values).reshape(max_number_of_nodes+1, max_number_of_nodes+1)


def join_nodes(n, mat, nodes_head):
    num_nodes = n
    num_active = n
    num_clusters = n

    joins = deque()

    while num_active > 3:
        if num_active == 4 and num_clusters == 2:
            p = nodes_head.next
            q = p.next if (p.next != p.nbr) else p.next.next
            if mat[p.id][q.id] + mat[p.nbr.id][q.nbr.id] < mat[p.id][q.nbr.id] + mat[p.nbr.id][q.id]:
                join3way(p, q, q.nbr, joins, mat, nodes_head, num_nodes)
            else:
                join3way(p, q.nbr, q, joins, mat, nodes_head, num_nodes)
            num_nodes += 2
            break

        p = nodes_head.next
        while p is not None:
            p.Sx = 0.0
            p = p.next

        p = nodes_head.next
        while p is not None:
            if p.nbr is None or p.nbr.id > p.id:
                q = p.next
                while q is not None:
                    if q.nbr is None or (q.nbr.id > q.id) and (q.nbr != p):
                        if p.nbr is None and q.nbr is None:
                            d_pq = mat[p.id][q.id]
                        elif p.nbr is not None and q.nbr is None:
                            d_pq = (mat[p.id][q.id] + mat[p.nbr.id][q.id]) / 2.0
                        elif p.nbr is None and q.nbr is not None:
                            d_pq = (mat[p.id][q.id] + mat[p.id][q.nbr.id]) / 2.0
                        else:
                            d_pq = (mat[p.id][q.id] + mat[p.id][q.nbr.id] + mat[p.nbr.id][q.id] + mat[p.nbr.id][
                                q.nbr.id]) / 4.0
                        p.Sx += d_pq
                        if p.nbr is not None:
                            p.nbr.Sx += d_pq
                        q.Sx += d_pq
                        if q.nbr is not None:
                            q.nbr.Sx += d_pq
                    q = q.next
            p = p.next

        Cx = None
        Cy = None
        best = 0.0

        p = nodes_head.next
        while p is not None:
            if p.nbr is not None and (p.nbr.id < p.id):
                p = p.next
                continue
            q = nodes_head.next
            while q is not None:
                if q == p:
                    break
                if (q.nbr is not None) and (q.nbr.id < q.id):
                    q = q.next
                    continue
                if q.nbr == p:
                    q = q.next
                    continue
                if (p.nbr is None) and (q.nbr is None):
                    d_pq = mat[p.id][q.id]
                elif (p.nbr is not None) and (q.nbr is None):
                    d_pq = (mat[p.id][q.id] + mat[p.nbr.id][q.id]) / 2.0
                elif (p.nbr is None) and (q.nbr is not None):
                    d_pq = (mat[p.id][q.id] + mat[p.id][q.nbr.id]) / 2.0
                else:
                    d_pq = (mat[p.id][q.id] + mat[p.id][q.nbr.id] + mat[p.nbr.id][q.id] + mat[p.nbr.id][q.nbr.id]) / 4.0

                q_pq = (num_clusters - 2.0) * d_pq - p.Sx - q.Sx

                if (Cx is None or (q_pq < best)) and (p.nbr != q):
                    Cx = p
                    Cy = q
                    best = q_pq
                q = q.next
            p = p.next

        x = Cx
        y = Cy

        if (Cx.nbr is not None) or (Cy.nbr is not None):
            Cx.Rx = compute_rx(Cx, Cx, Cy, mat, nodes_head)
            if Cx.nbr is not None:
                Cx.nbr.Rx = compute_rx(Cx.nbr, Cx, Cy, mat, nodes_head)
            Cy.Rx = compute_rx(Cy, Cx, Cy, mat, nodes_head)
            if Cy.nbr is not None:
                Cy.nbr.Rx = compute_rx(Cy.nbr, Cx, Cy, mat, nodes_head)

        m = num_clusters
        if Cx.nbr is not None:
            m += 1
        if Cy.nbr is not None:
            m += 1

        best = (m - 2.0) * mat[Cx.id][Cy.id] - Cx.Rx - Cy.Rx
        if Cx.nbr is not None:
            q_pq = (m - 2.0) * mat[Cx.nbr.id][Cy.id] - Cx.nbr.Rx - Cy.Rx
            if q_pq < best:
                x = Cx.nbr
                y = Cy
                best = q_pq

        if Cy.nbr is not None:
            q_pq = (m - 2.0) * mat[Cx.id][Cy.nbr.id] - Cx.Rx - Cy.nbr.Rx
            if q_pq < best:
                x = Cx
                y = Cy.nbr
                best = q_pq

        if (Cx.nbr is not None) and (Cy.nbr is not None):
            q_pq = (m - 2.0) * mat[Cx.nbr.id][Cy.nbr.id] - Cx.nbr.Rx - Cy.nbr.Rx
            if q_pq < best:
                x = Cx.nbr
                y = Cy.nbr

        if x.nbr is None and y.nbr is None:
            join2way(x, y)
            num_clusters -= 1
        elif x.nbr is None:
            join3way(x, y, y.nbr, joins, mat, nodes_head, num_nodes)
            num_nodes += 2
            num_active -= 1
            num_clusters -= 1
        elif (y.nbr is None) or (num_active == 4):
            join3way(y, x, x.nbr, joins, mat, nodes_head, num_nodes)
            num_nodes += 2
            num_active -= 1
            num_clusters -= 1
        else:
            num_nodes = join4way(x.nbr, x, y, y.nbr, joins, mat, nodes_head, num_nodes)
            num_active -= 2
            num_clusters -= 1

    return joins


def join2way(x, y):
    x.nbr = y
    y.nbr = x


def join3way(x, y, z, joins, mat, nodes_head, num_nodes):
    u = NetNode(num_nodes + 1)
    u.ch1 = x
    u.ch2 = y

    v = NetNode(num_nodes + 2)
    v.ch1 = y
    v.ch2 = z

    u.next = x.next
    u.prev = x.prev

    if u.next is not None:
        u.next.prev = u
    if u.prev is not None:
        u.prev.next = u

    v.next = z.next
    v.prev = z.prev
    if v.next is not None:
        v.next.prev = v
    if v.prev is not None:
        v.prev.next = v

    if y.next is not None:
        y.next.prev = y.prev
    if y.prev is not None:
        y.prev.next = y.next

    u.nbr = v
    v.nbr = u

    p = nodes_head.next
    while p is not None:
        mat[u.id][p.id] = mat[p.id][u.id] = (2.0 / 3.0) * mat[x.id][p.id] + mat[y.id][p.id] / 3.0
        mat[v.id][p.id] = mat[p.id][v.id] = (2.0 / 3.0) * mat[z.id][p.id] + mat[y.id][p.id] / 3.0
        p = p.next
    mat[u.id][u.id] = mat[v.id][v.id] = 0.0

    joins.append(u)

    return u


def join4way(x2, x, y, y2, joins, mat, nodes_head, num_nodes):
    u = join3way(x2, x, y, joins, mat, nodes_head, num_nodes)
    num_nodes += 2
    join3way(u, u.nbr, y2, joins, mat, nodes_head, num_nodes)
    num_nodes += 2
    return num_nodes


def compute_rx(z, Cx, Cy, mat, nodes_head):
    Rx = 0.0

    p = nodes_head.next
    while p is not None:
        if p == Cx or p == Cx.nbr or p == Cy or p == Cy.nbr or p.nbr is None:
            Rx += mat[z.id][p.id]
        else:
            Rx += mat[z.id][p.id] / 2.0
        p = p.next
    return Rx


def expand_nodes(joins, nodes_head):
    x = nodes_head.next
    y = x.next
    z = y.next
    z.next = x
    x.prev = z

    while len(joins) > 0:
        u = joins.pop()
        v = u.nbr
        x = u.ch1
        y = u.ch2
        z = v.ch2
        if v != u.next:
            tmp = u
            u = v
            v = tmp
            tmp = x
            x = z
            z = tmp

        x.prev = u.prev
        x.prev.next = x
        x.next = y
        y.prev = x
        y.next = z
        z.prev = y
        z.next = v.next
        z.next.prev = z

    while x.id != 1:
        x = x.next

    cycle = []
    a = x
    while True:
        cycle.append(a.id)
        a = a.next
        if a == x:
            break

    return cycle


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

    def create_string(self):
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



