# nnet_splits.py
"""Runs the neighbor-net algorithm to compute circular splits

See: Bryant and Moulton (2004)
See: Huson and Bryant (2006)


LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""

from splitspy.splits.basic_split import *
import math

__author__ = "David J. Bryant and Daniel H. Huson"

CG_EPSILON = 0.0001


def compute(n_tax: int, mat: [float], cycle: [int], cutoff=0.00001, constrained=True) -> [Split]:
    if n_tax == 1:
        return []
    elif n_tax == 2:
        return [cyc_split([1, 2], 2, 2, mat[1][2])] if mat[1][2] >= cutoff else []

    d = setup_d(n_tax, mat, cycle)
    x = []

    if not constrained:
        unconstrained_least_squares(n_tax, d, x)
    else:
        w = setup_w(n_tax)
        x = active_conjugate(n_tax, d, w)

    splits = []

    index = 0
    for i in range(1, n_tax + 1):
        for j in range(i + 1, n_tax + 1):
            if x[index] > cutoff:
                splits.append(cyc_split(cycle, i + 1, j, x[index]))
            index += 1

    return splits


def setup_d(n: int, mat: [float], cycle: [int]) -> [float]:
    d = []

    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            d.append(mat[cycle[i] - 1][cycle[j] - 1])
    return d


def setup_w(n: int) -> [float]:
    w = []
    for i in range(0, n):
        for j in range(i + 1, n):
            w.append(1.0)
    return w


def unconstrained_least_squares(n_tax: int, d: [float], x: [float]) -> None:
    x.clear()
    index = 0
    for i in range(0, n_tax - 3 + 1):
        x.append((d[index] + d[index + (n_tax - i - 2) + 1] - d[index + 1]) / 2.0)
        index += 1
        for j in range(i + 2, n_tax - 2 + 1):
            x.append((d[index] + d[index + (n_tax - i - 2) + 1] - d[index + 1] - d[index + (n_tax - i - 2)]) / 2.0)
            index += 1
        if i == 0:
            x.append((d[0] + d[n_tax - 2] - d[2 * n_tax - 4]) / 2.0)
        else:
            x.append((d[index] + d[i] - d[i - 1] - d[index + (n_tax - i - 2)]) / 2.0)
        index += 1
    x.append((d[index] + d[n_tax - 2] - d[n_tax - 3]) / 2.0)


def active_conjugate(n_tax: int, d: [float], w: [float]) -> [float]:
    x: [float]
    x = []

    unconstrained_least_squares(n_tax, d, x)

    if all(a >= 0 for a in x):
        return x

    n_pairs = len(d)

    active = [False] * n_pairs

    y = [0.0] * n_pairs
    for k in range(0, n_pairs):
        y[k] = w[k] * d[k]

    at_wd = calculate_Atx(n_tax, y)

    old_x = [1.0] * n_pairs

    first_pass = True

    while True:
        while True:
            if first_pass:
                first_pass = False
            else:
                circular_conjugate_grads(n_tax, n_pairs, w, at_wd, active, x)

            to_contract = worst_indices(x, 0.6)
            if len(to_contract) > 0:
                for index in to_contract:
                    x[index] = 0.0
                    active[index] = True
                circular_conjugate_grads(n_tax, n_pairs, w, at_wd, active, x)

            min_i = -1
            min_xi = -1.0
            for i in range(0, n_pairs):
                if x[i] < 0.0:
                    xi = (old_x[i]) / (old_x[i] - x[i])
                    if min_i == -1 or xi < min_xi:
                        min_i = i
                        min_xi = xi
            if min_i == -1:
                break
            else:
                for i in range(0, n_pairs):
                    if not active[i]:
                        old_x[i] += min_xi * (x[i] - old_x[i])
                active[min_i] = True
                x[min_i] = 0.0

        y = calculate_AB(n_tax, x)
        for i in range(0, n_pairs):
            y[i] *= w[i]
        r = calculate_Atx(n_tax, y)

        min_i = -1
        min_grad = 1.0

        for i in range(0, n_pairs):
            r[i] -= at_wd[i]
            r[i] *= 2.0
            if active[i]:
                grad_ij = r[i]
                if min_i == -1 or grad_ij < min_grad:
                    min_i = i
                    min_grad = grad_ij

        if min_i == -1 or min_grad > -0.0001:
            break
        else:
            active[min_i] = False
    return x


def calculate_Atx(n: int, d: [float]) -> [float]:
    p = [0.0] * len(d)

    index = 0
    for i in range(0, n - 1):
        p[index] = row_sum(n, d, i + 1)
        index += (n - i - 1)

    index = 1
    for i in range(0, n - 2):
        p[index] = p[index - 1] + p[index + (n - i - 2)] - 2 * d[index + (n - i - 2)]
        index += (n - i - 2) + 1

    for k in range(3, n):
        index = k - 1
        for i in range(0, n - k):
            p[index] = p[index - 1] + p[index + n - i - 2] - p[index + n - i - 3] - 2.0 * d[index + n - i - 2]
            index += (n - i - 2) + 1
    return p


def row_sum(n: int, d: [float], k: int) -> float:
    r = 0
    index = 0

    if k > 0:
        index = k - 1
        for i in range(0, k):
            r += d[index]
            index += (n - i - 2)
        index += 1

    for j in range(k + 1, n):
        r += d[index]
        index += 1

    return r


def worst_indices(x: [float], prop_kept: float) -> [int]:
    if prop_kept == 0.0:
        return []

    prop_kept = 0.1

    n_pairs = len(x)

    x_cpy = []
    for value in x:
        if value < 0:
            x_cpy.append(value)
    n_neg = len(x_cpy)

    if n_neg == 0:
        return []

    x_cpy.sort()

    n_kept = math.ceil(prop_kept * n_neg)
    cutoff = x_cpy[n_kept - 1]

    front = 0
    back = n_kept - 1

    worst = [0] * n_kept

    for i in range(0, n_pairs):
        if x[i] < cutoff:
            worst[front] = i
            front += 1
        elif x[i] == cutoff:
            if back >= front:
                worst[back] = i
                back -= 1

    return worst


def circular_conjugate_grads(n_tax: int, n_pairs: int, W: [float], b: [float], active: [bool], x: [float]) -> None:
    k_max = n_tax * (n_tax - 1) / 2

    y = calculate_AB(n_tax, x)

    for k in range(0, n_pairs):
        y[k] = W[k] * y[k]

    r = calculate_Atx(n_tax, y)

    for k in range(0, n_pairs):
        if not active[k]:
            r[k] = b[k] - r[k]
        else:
            r[k] = 0.0

    rho = norm(r)
    rho_old = 0

    e_0 = CG_EPSILON * math.sqrt(norm(b))
    k = 0

    p = []
    while rho > e_0 * e_0 and k < k_max:
        k = k + 1
        if k == 1:
            p = r.copy()
        else:
            beta = rho / rho_old
            for i in range(0, n_pairs):
                p[i] = r[i] + beta * p[i]

        y = calculate_AB(n_tax, p)

        for i in range(0, n_pairs):
            y[i] *= W[i]

        u = calculate_Atx(n_tax, y)

        for i in range(0, n_pairs):
            if active[i]:
                u[i] = 0.0

        alpha = 0.0
        for i in range(0, n_pairs):
            alpha += p[i] * u[i]

        alpha = rho / alpha

        for i in range(0, n_pairs):
            x[i] += alpha * p[i]
            r[i] -= alpha * u[i]

        rho_old = rho
        rho = norm(r)


def calculate_AB(n: int, b: [float]) -> [float]:
    d = [0.0] * len(b)

    d_index = 0
    for i in range(0, n - 1):
        d_ij = 0.0
        index = i - 1
        for k in range(0, i):
            d_ij += b[index]
            index += n - k - 2
        index += 1
        for k in range(i + 1, n):
            d_ij += b[index]
            index += 1

        d[d_index] = d_ij
        d_index += (n - i - 2) + 1

    index = 1
    for i in range(0, n - 2):
        d[index] = d[index - 1] + d[index + (n - i - 2)] - 2 * b[index - 1]
        index += 1 + (n - i - 2)

    for k in range(3, n):
        index = k - 1
        for i in range(0, n - k):
            d[index] = d[index - 1] + d[index + (n - i - 2)] - d[index + (n - i - 2) - 1] - 2.0 * b[index - 1]
            index += 1 + (n - i - 2)

    return d


def norm(x: [float]) -> float:
    n = 0.0
    for value in x:
        n += value * value
    return n
