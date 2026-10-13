"""Enriched test subspace for cubic graphs with a cut vertex.

W(m, k) = span{ 1_{C_i} - n_i e_z : i = 1,2,3 }  (+)  span{ e_u - (1/m)1_{C_m} :
u in C_m }  (+)  span{ e_u - (1/k)1_{C_k} : u in C_k }.

Every vector in W is orthogonal to the all-ones vector, so
    lambda_2(A(G))  >=  max_{x in W, x != 0}  x'A(G)x / ||x||^2 .
The bound depends only on the arm orders (m, k, n_3) via this construction.
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


def build(m, k, l):
    off = 0
    edges, tips, orders = [], [], []
    for o in (m, k, l):
        n, e, t = arm_of_order(o)
        edges += [(u + off, v + off) for u, v in e]
        tips.append(off + t)
        orders.append(n)
        off += n
    z = off
    edges += [(t, z) for t in tips]
    n = z + 1
    A = np.zeros((n, n))
    for u, v in edges:
        A[u, v] = A[v, u] = 1
    return n, A, orders, [0, m, m + k], z


def enriched_bound(m, k, l, exact_check=False):
    n, A, orders, offs, z = build(m, k, l)
    cols = []
    # arm indicators
    for i in range(3):
        v = np.zeros(n)
        v[offs[i]:offs[i] + orders[i]] = 1
        v[z] = -orders[i]
        cols.append(v)
    # mean-zero directions inside arms 0 and 1
    for i in (0, 1):
        for u in range(orders[i]):
            v = np.zeros(n)
            v[offs[i] + u] += 1
            v[offs[i]:offs[i] + orders[i]] -= 1.0 / orders[i]
            cols.append(v)
    V = np.column_stack(cols)
    U, s, _ = np.linalg.svd(V, full_matrices=False)
    r = int((s > 1e-9).sum())
    Q = U[:, :r]
    assert abs(Q.sum(axis=0)).max() < 1e-8, "Q must be orthogonal to 1"
    B = Q.T @ (A @ Q)
    val = float(np.linalg.eigvalsh((B + B.T) / 2).max())
    if exact_check:
        lam2 = float(np.linalg.eigvalsh(A)[::-1][1])
        return val, lam2
    return val


def main():
    print(f"2*sqrt(2) = {RHO3:.12f}\n")
    print(f"{'arms':>20}{'bound':>14}{'exact lam2':>14}   ok?   >2sqrt2?")
    for (m, k, l) in [(5, 5, 5), (5, 5, 7), (5, 5, 9), (5, 5, 21), (5, 5, 101),
                      (5, 5, 1001), (5, 7, 5), (5, 7, 101), (7, 7, 7), (9, 9, 9),
                      (5, 23, 5), (5, 23, 101), (21, 23, 23), (5, 199, 199),
                      (9, 11, 13), (11, 13, 17)]:
        b, lam2 = enriched_bound(m, k, l, exact_check=True)
        ok = b <= lam2 + 1e-9
        print(f"{str((m, k, l)):>20}{b:14.9f}{lam2:14.9f}   {str(ok):>5}"
              f"   {'YES' if b > RHO3 else 'NO'}")


if __name__ == "__main__":
    main()