"""Closed form of the 6-dimensional order-dependent pencil, and verification.

For a connected cubic graph G with cut vertex z and G - z = C_1 U C_2 U C_3,
each arm C_i of odd order n_i >= 5 meeting z at a single vertex t_i, put

    p_i = 1_{C_i} - n_i e_z,      q_i = e_{t_i} - (1/n_i) 1_{C_i}.

Both are orthogonal to the all-ones vector.  Writing
V = [p_1 q_1 p_2 q_2 p_3 q_3] and M = V'AV, K = V'V, one gets (i, j index arms):

    M[p_i,p_i] = n_i - 1,                M[p_i,p_j] = 0            (i != j)
    M[p_i,q_i] = -(n_i-1) - (n_i-1)/n_i, M[p_i,q_j] = -n_i (n_j-1)/n_j
    M[q_i,q_i] = (2n_i^2-n_i-1)/n_i^2,   M[q_i,q_j] = 0
    K[p_i,p_i] = n_i + n_i^2,            K[p_i,p_j] = n_i n_j - n_i - n_j
    K[p_i,q_j] = 0,                      K[q_i,q_i] = (n_i-1)^2/n_i^2

lambda_2(A(G)) >= lambda_max(M, K).
"""
from __future__ import annotations

import itertools

import numpy as np

RHO3 = 2 * np.sqrt(2)


def pencil(n):
    M = np.zeros((6, 6))
    K = np.zeros((6, 6))
    idx = {i: (2 * i, 2 * i + 1) for i in range(3)}
    for i, ni in enumerate(n):
        for j, nj in enumerate(n):
            pi, qi = idx[i]
            pj, qj = idx[j]
            M[pi, pj] = (ni - 1) if i == j else 0.0
            M[qi, qj] = ((2 * ni * ni - ni - 1) / ni ** 2) if i == j else 0.0
            M[pi, qj] = (-(ni - 1) - (ni - 1) / ni) if i == j \
                else -ni * (nj - 1) / nj
            M[pj, qi] = M[pi, qj]
            K[pi, pj] = (ni + ni * ni) if i == j else (ni * nj - ni - nj)
            K[qi, qj] = ((ni - 1) / ni) ** 2 if i == j else 0.0
    return (M + M.T) / 2, (K + K.T) / 2


def gen_pencil(n):
    M, K = pencil(n)
    ev, U = np.linalg.eigh(K)
    Ki = U @ np.diag(1.0 / np.sqrt(ev)) @ U.T
    S = Ki @ M @ Ki
    S = (S + S.T) / 2
    w, Vv = np.linalg.eigh(S)
    a = Ki @ Vv[:, -1]
    return float(w[-1]), a


# ---- explicit-vector verification -------------------------------------------------
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


def explicit(n_orders):
    off = 0
    edges, tips = [], []
    for o in n_orders:
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
    for o in n_orders:
        offs.append(o2)
        o2 += o
    V = []
    for i, o in enumerate(n_orders):
        p = np.zeros(N)
        p[offs[i]:offs[i] + o] = 1
        p[z] = -o
        V.append(p)
        q = np.zeros(N)
        q[tips[i]] += 1
        q[offs[i]:offs[i] + o] -= 1.0 / o
        V.append(q)
    return A, np.column_stack(V)


def explicit_bound(n_orders):
    A, V = explicit(n_orders)
    K = V.T @ V
    M = V.T @ (A @ V)
    ev, U = np.linalg.eigh(K)
    Ki = U @ np.diag(1.0 / np.sqrt(ev)) @ U.T
    S = Ki @ M @ Ki
    return float(np.linalg.eigvalsh((S + S.T) / 2).max())


def main():
    print("cross-check: closed-form pencil vs explicit coordinate vectors")
    for n in [(5, 5, 5), (5, 7, 9), (5, 5, 11), (7, 9, 11), (5, 23, 5)]:
        b1, _ = gen_pencil(n)
        b2 = explicit_bound(n)
        print(f"  {str(n):>12}: closed form {b1:.12f}  explicit {b2:.12f}  "
              f"{'OK' if abs(b1-b2) < 1e-9 else 'MISMATCH'}")

    print(f"\n2*sqrt(2) = {RHO3:.12f}")
    print("\nsmallest bound over each family (arm orders):")
    for n in [(5, 5, 5), (5, 5, 7), (5, 5, 11), (5, 5, 51), (5, 5, 1001),
              (5, 7, 7), (5, 9, 9), (7, 7, 7), (5, 23, 5), (9, 11, 13),
              (11, 13, 17), (5, 999, 999), (999, 1001, 1003)]:
        b, _ = gen_pencil(n)
        print(f"  {str(n):>18}: bound = {b:.9f}  {'YES' if b > RHO3 else 'NO'}")


if __name__ == "__main__":
    main()