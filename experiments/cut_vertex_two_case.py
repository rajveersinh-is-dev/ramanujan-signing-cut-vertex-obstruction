"""The two-case exclusion argument for cubic graphs with a cut vertex.

Setup.  G connected cubic, z a cut vertex, G - z = C_1 U C_2 U C_3, arm orders
odd, n_1 <= n_2 <= n_3, each arm joined to z at a single vertex t_i.
Define p_i = 1_{C_i} - n_i e_z  and  q_i = e_{t_i} - (1/n_i) 1_{C_i}; all are
orthogonal to the all-ones vector.

Case A (n_2 >= 7):  the two-arm indicator test
    x = (1/n_i) 1_{C_i} - (1/n_j) 1_{C_j}  gives
    lambda_2 >= 3 - (n_i^2+n_j^2)/(n_i n_j (n_i+n_j))
             = 3 - psi(r)/n_i,   psi(r) = (r+1/r)/(1+r) <= 1 for r >= 1,
    so lambda_2 >= 3 - 1/n_2 >= 20/7 > 2 sqrt 2.

Case B (n_2 = 5, hence n_1 = 5):  the 2-dimensional subspace
    span{ p_1 + p_2 , q_1 + q_2 }  gives a bound depending only on n_1 = n_2 = 5.

Everything below is computed from explicit coordinate vectors on a concrete
cubic graph, and cross-checked against exact spectra.
"""
from __future__ import annotations

import numpy as np

RHO3 = 2 * np.sqrt(2)


def arm_of_order(n):
    """Concrete arm of odd order n: cubic except at the tip (degree 2)."""
    assert n % 2 == 1 and n >= 5
    if n == 5:
        he = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    else:
        kk = (n - 1) // 2
        he = sorted(set([(i, (i + 1) % kk) for i in range(kk)]
                        + [(i + kk, ((i + 1) % kk) + kk) for i in range(kk)]
                        + [(i, i + kk) for i in range(kk)]))
    E = [e for e in he if frozenset(e) != frozenset((0, 1))]
    return n, sorted(E + [(n - 1, 0), (n - 1, 1)]), n - 1


def build(orders):
    off, edges, tips = 0, [], []
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
    return N, A, offs, z, tips


def p_vec(N, offs, z, i, n):
    v = np.zeros(N)
    v[offs[i]:offs[i] + n] = 1
    v[z] = -n
    return v


def q_vec(N, offs, z, i, n, tip):
    v = np.zeros(N)
    v[tip] += 1
    v[offs[i]:offs[i] + n] -= 1.0 / n
    return v


def top_of_span(A, vectors):
    V = np.column_stack(vectors)
    assert abs(V.sum(axis=0)).max() < 1e-8
    K = V.T @ V
    M = V.T @ (A @ V)
    ev, U = np.linalg.eigh(K)
    Ki = U @ np.diag(1.0 / np.sqrt(ev)) @ U.T
    S = (Ki @ M @ Ki)
    return float(np.linalg.eigvalsh((S + S.T) / 2).max())


def case_b_bound(n1, n2):
    """Bound from span{p_1+p_2, q_1+q_2} for a graph whose two arms have
    orders n1 = n2 = 5 and third arm of order l (any)."""
    l = 5
    N, A, offs, z, tips = build((n1, n2, l))
    u = p_vec(N, offs, z, 0, n1) + p_vec(N, offs, z, 1, n2)
    v = q_vec(N, offs, z, 0, n1, tips[0]) + q_vec(N, offs, z, 1, n2, tips[1])
    return top_of_span(A, [u, v])


def case_a_bound(n2, n3):
    """Bound from the two-arm indicator test on arms 2 and 3."""
    N, A, offs, z, tips = build((5, n2, n3))
    x = (p_vec(N, offs, z, 1, n2) + (n2 * 0) * 0)  # placeholder
    # the test vector is  (1/n2) 1_{C2} - (1/n3) 1_{C3}
    x = np.zeros(N)
    x[offs[1]:offs[1] + n2] = 1.0 / n2
    x[offs[2]:offs[2] + n3] = -1.0 / n3
    return float(x @ A @ x / (x @ x))


def psi(r):
    return (r + 1.0 / r) / (1.0 + r)


def main():
    print(f"2*sqrt(2) = {RHO3:.12f}")
    print(f"20/7      = {20/7:.12f}")

    print("\nCase B: span{ p_1+p_2 , q_1+q_2 } with n_1 = n_2 = 5")
    b = case_b_bound(5, 5)
    print(f"  bound = {b:.12f}   > 2 sqrt 2 ? {b > RHO3}")
    print(f"  (independent of the third arm: verify)")
    for l in (5, 7, 9, 101, 1001):
        N, A, offs, z, tips = build((5, 5, l))
        u = p_vec(N, offs, z, 0, 5) + p_vec(N, offs, z, 1, 5)
        v = q_vec(N, offs, z, 0, 5, tips[0]) + q_vec(N, offs, z, 1, 5, tips[1])
        bb = top_of_span(A, [u, v])
        lam2 = float(np.linalg.eigvalsh(A)[::-1][1])
        print(f"    l={l:5d}: bound={bb:.12f}  lam2(graph)={lam2:.9f}  "
              f"{'OK' if bb <= lam2 + 1e-9 else 'VIOLATION'}")

    print("\nCase A: two-arm indicator test on arms 2,3 (n_2 >= 7)")
    for (n2, n3) in [(7, 7), (7, 9), (7, 101), (9, 9), (21, 23), (199, 201)]:
        a = case_a_bound(n2, n3)
        N, A, offs, z, tips = build((5, n2, n3))
        lam2 = float(np.linalg.eigvalsh(A)[::-1][1])
        print(f"  n2={n2:4d} n3={n3:4d}: bound={a:.9f}  lam2={lam2:.9f}  "
              f"{'OK' if a <= lam2 + 1e-9 else 'VIOLATION'}")
    print(f"  formula 3 - psi(n3/n2)/n2 >= 3 - 1/n2 >= 20/7 = {20/7:.9f}")
    print(f"  max of psi over r>=1 is 1 (attained at r=1): "
          f"psi(1) = {psi(1.0):.6f}")


if __name__ == "__main__":
    main()