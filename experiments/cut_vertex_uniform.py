"""Uniform exclusion test for connected cubic graphs with a cut vertex.

Setup.  G connected cubic, z a cut vertex, G - z = C_1 U C_2 U C_3, each arm
C_i of odd order n_i >= 5 joined to z at a single vertex t_i.

Test subspace.  For two arms with orders m <= k (m, k odd, >= 5) consider
x = a*(1_{C_m} - m*e_z) + b*(1_{C_k} - k*e_z) + xi*(e_{t_m} - (1/m) 1_{C_m})
    + eta*(e_{t_k} - (1/k) 1_{C_k}),
supported on {z} U C_m U C_k and vanishing on the third arm.  Every such x
satisfies <x, 1> = 0, hence lambda_2(A(G)) >= max_{x != 0} x'A(G)x / ||x||^2
over this 4-dimensional subspace W.

Everything below is computed from explicit coordinate vectors: no hand-derived
matrix is used, and the answer is cross-checked against the exact lambda_2 of
the concrete cubic graph.
"""
from __future__ import annotations

import numpy as np

RHO3 = 2 * np.sqrt(2)


def arm_of_order(n):
    """Concrete arm of odd order n: cubic except one degree-2 vertex (the tip)."""
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


def sample_graph(m, k, l=5):
    """Cubic graph: three arms of orders m, k, l joined at a centre z."""
    off = 0
    edges, tips = [], []
    for o in (m, k, l):
        n, e, t = arm_of_order(o)
        edges += [(u + off, v + off) for u, v in e]
        tips.append(off + t)
        off += n
    z = off
    edges += [(t, z) for t in tips]
    return z + 1, sorted(edges)


def test_subspace(m, k, l=5):
    """Basis (list of n-vectors) of the 4-dimensional test subspace W."""
    n, edges = sample_graph(m, k, l)
    A = np.zeros((n, n))
    for u, v in edges:
        A[u, v] = A[v, u] = 1
    off_m = 0
    off_k = m
    z = m + k + l
    tm = 0            # tip of first arm is vertex 0 of its arm
    tk = off_k + 0
    e = np.zeros(n)
    v1 = np.zeros(n); v1[off_m:off_m + m] = 1; v1[z] = -m
    v2 = np.zeros(n); v2[off_k:off_k + k] = 1; v2[z] = -k
    v3 = np.zeros(n); v3[tm] += 1; v3[off_m:off_m + m] -= 1.0 / m
    v4 = np.zeros(n); v4[tk] += 1; v4[off_k:off_k + k] -= 1.0 / k
    V = np.column_stack([v1, v2, v3, v4])
    assert np.allclose(V.sum(axis=0), 0), V.sum(axis=0)
    return n, edges, A, V


def subspace_bound(m, k, l=5):
    n, edges, A, V = test_subspace(m, k, l)
    K = V.T @ V
    M = V.T @ (A @ V)
    ev, U = np.linalg.eigh(K)
    Ki = U @ np.diag(1.0 / np.sqrt(ev)) @ U.T
    S = Ki @ M @ Ki
    return float(np.linalg.eigvalsh((S + S.T) / 2).max())


def bottleneck(m, k):
    """Two-arm indicator test: 3 - (m^2+k^2)/(mk(m+k))."""
    return 3.0 - (m * m + k * k) / float(m * k * (m + k))


def main():
    print("2*sqrt(2) =", f"{RHO3:.12f}")
    print("\ncross-check: subspace bound must not exceed the exact lambda_2 of "
          "a concrete graph with the same arm orders")
    for (m, k) in [(5, 5), (5, 7), (5, 9), (5, 23), (7, 7), (9, 9), (21, 23)]:
        b = subspace_bound(m, k)
        n, edges = sample_graph(m, k, max(m, k))
        lam2 = float(np.linalg.eigvalsh(
            np.array([[1.0 if (u, v) in set(edges) or (v, u) in set(edges)
                       else 0.0 for v in range(n)] for u in range(n)]))[::-1][1])
        ok = b <= lam2 + 1e-9
        print(f"  m={m:3d} k={k:3d}  bound={b:.9f}  lam2(sample)={lam2:.9f}"
              f"  {'OK' if ok else 'VIOLATION'}")

    print("\nscan: min of the subspace bound over odd 5 <= m <= k <= 600")
    odds = list(range(5, 601, 2))
    worst, at = 1e18, None
    bad = []
    for i, m in enumerate(odds):
        for k in odds[i:]:
            b = subspace_bound(m, k)
            if b < worst:
                worst, at = b, (m, k)
            if b <= RHO3 + 1e-12:
                bad.append((m, k, b))
    print(f"  min = {worst:.12f} at (m,k) = {at}")
    print(f"  #pairs with bound <= 2sqrt2: {len(bad)}")
    for x in bad[:20]:
        print("   ", x)
    print()
    for (m, k) in [(5, 5), (5, 7), (5, 9), (5, 23), (5, 51), (5, 199), (7, 7),
                   (7, 9), (9, 9), (21, 23), (199, 201), (599, 601)]:
        b = subspace_bound(m, k)
        print(f"  m={m:4d} k={k:4d}  bottleneck={bottleneck(m,k):.9f}  "
              f"bound={b:.9f}  {'YES' if b > RHO3 else 'NO'}")


if __name__ == "__main__":
    main()