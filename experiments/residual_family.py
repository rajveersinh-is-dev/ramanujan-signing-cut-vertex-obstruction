"""The residual family: arm orders (5, 5, l).

The two-arm indicator (bottleneck) test gives
    lambda_2 >= 3 - (m^2+k^2)/(mk(m+k)),
which for m = 5 exceeds 2 sqrt 2 exactly when  k/m in (rho_-, rho_+) with
rho_\pm the roots of psi(r) = 5 (3 - 2 sqrt 2), psi(r) = (r + 1/r)/(1 + r);
in particular it exceeds 2 sqrt 2 whenever the smaller arm order is >= 7.
So the only residual case is two arms of order 5, handled here with the
6-dimensional subspace

    W = span{ 1_{C_i} - n_i e_z ,  e_{t_i} - (1/n_i) 1_{C_i} : i = 1,2,3 }.

This script prints the top Rayleigh quotient and the top eigenvector of the
associated 6x6 pencil, looking for a closed-form test vector.
"""
from __future__ import annotations

import numpy as np

RHO3 = 2 * np.sqrt(2)


def arm_of_order(n):
    assert n % 2 == 1 and n >= 5
    if n == 5:
        he = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    else:
        kk = (n - 1) // 2
        he = sorted(set([(i, (i + 1) % kk) for i in range(kk)]
                        + [(i + kk, ((i + 1) % kk) + kk) for i in range(kk)]
                        + [(i, i + kk) for i in range(kk)]))
    E = [e for e in he if frozenset(e) != frozenset((0, 1))]
    E = E + [(n - 1, 0), (n - 1, 1)]
    return n, sorted(E), n - 1


def pencil_data(orders):
    off = 0
    edges, tips = [], []
    for o in orders:
        nn, e, t = arm_of_order(o)
        edges += [(u + off, v + off) for u, v in e]
        tips.append(off + t)
        off += nn
    z = off
    edges += [(t, z) for t in tips]
    N = z + 1
    A = np.zeros((N, N))
    for u, v in edges:
        A[u, v] = A[v, u] = 1
    offs, o2 = [], 0
    for o in orders:
        offs.append(o2)
        o2 += o
    V = []
    for i, o in enumerate(orders):
        p = np.zeros(N)
        p[offs[i]:offs[i] + o] = 1
        p[z] = -o
        V.append(p)
        q = np.zeros(N)
        q[tips[i]] += 1
        q[offs[i]:offs[i] + o] -= 1.0 / o
        V.append(q)
    V = np.column_stack(V)
    K = V.T @ V
    M = V.T @ (A @ V)
    ev, U = np.linalg.eigh(K)
    Ki = U @ np.diag(1.0 / np.sqrt(ev)) @ U.T
    S = (Ki @ M @ Ki)
    S = (S + S.T) / 2
    w, Vv = np.linalg.eigh(S)
    return float(w[-1]), Ki @ Vv[:, -1]


def psi(r):
    return (r + 1.0 / r) / (1.0 + r)


def rho_pm(c):
    # psi(r) = c  ->  (1+r^2) = c r (1+r)  ->  (1-c) r^2 - c r + 1 = 0
    a, b, cc = (1 - c), (-c), 1.0
    disc = b * b - 4 * a * cc
    return ((-b - np.sqrt(disc)) / (2 * a), (-b + np.sqrt(disc)) / (2 * a))


def main():
    c = 5 * (3 - 2 * np.sqrt(2))
    lo, hi = rho_pm(c)
    print(f"psi(r) = 5(3-2sqrt2) = {c:.9f}  ->  r in ({lo:.6f}, {hi:.6f})")
    print(f"  so the bottleneck test works for 5 < r < 4.4576..., i.e. "
          f"odd k in {{9,...,21}}\n")
    print(f"2*sqrt(2) = {RHO3:.12f}")
    print()
    print(f"{'l':>9}{'bound':>16}{'excess':>14}   top eigenvector "
          f"[a1,b1,c1,a2,b2,c2] (order p1,q1,p2,q2,p3,q3)")
    for l in (5, 7, 9, 11, 15, 21, 25, 31, 51, 101, 201, 501, 1001, 5001,
              10001, 20001):
        b, a = pencil_data((5, 5, l))
        s = " ".join(f"{x:8.4f}" for x in a)
        print(f"{l:9d}{b:16.9f}{b - RHO3:14.3e}   {s}")


if __name__ == "__main__":
    main()