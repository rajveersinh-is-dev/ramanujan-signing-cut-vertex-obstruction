"""Symbolic pencil for the residual case n_1 = n_2 = 5, from correct formulas.

Combinatorial facts used.  Each arm C_i has order n_i, all its vertices have
degree 3 except the tip t_i, which has degree 2 and is joined to z; hence t_i has
exactly two neighbours inside C_i.

Represent a vector by (z-value, [(arm, c, d), ...]) meaning: value c on
C_i \\ {t_i}, value c + d at t_i.  Then

  (1)  <x, y> |_C_i = (n_i - 1) c_1 c_2 + (c_1+d_1)(c_2+d_2)

  (2)  x' A y |_C_i = sum_u x_u sum_{v ~ u, v in C_i} y_v
                     = 3 (n_i-1) c_1 c_2 + 2 d_1 c_2 + 2 c_1 (c_2 + d_2)

       (group by v != t_i:  each contributes c_1 * d(v) y_v = 3 c_1 c_2; then
        the tip neighbour d_1 * y_v for v ~ t_i, v != t_i, giving 2 d_1 c_2;
        then y_{t_i} * sum_{u ~ t_i} x_u = y_{t_i} * 2 c_1.)

  (3)  the three edges z t_i contribute  sum_i [ z_1 (c_2+d_2) + (c_1+d_1) z_2 ]

The resulting pencil is cross-validated against explicitly built graphs.
"""
from __future__ import annotations

import importlib.util

import numpy as np
import sympy as sp

l = sp.symbols('l', positive=True)
S2 = sp.sqrt(2)
TAU = 2 * S2


def bilinear(g1, g2, sizes):
    (z1, a1), (z2, a2) = g1, g2
    d1 = {i: (c, dd) for i, c, dd in a1}
    d2 = {i: (c, dd) for i, c, dd in a2}
    gram = z1 * z2
    adj = sp.Integer(0)
    for i in sorted(sizes):
        n = sizes[i]
        c1, e1 = d1.get(i, (sp.Integer(0), sp.Integer(0)))
        c2, e2 = d2.get(i, (sp.Integer(0), sp.Integer(0)))
        t1, t2 = c1 + e1, c2 + e2
        gram += (n - 1) * c1 * c2 + t1 * t2
        adj += 3 * (n - 1) * c1 * c2 + 2 * e1 * c2 + 2 * c1 * t2
        adj += z1 * t2 + t1 * z2
    return sp.expand(gram), sp.expand(adj)


def pencil(n1, n2, n3):
    sizes = {1: n1, 2: n2, 3: n3}
    gens = [
        (sp.Integer(-(n1 + n2)),
         [(1, sp.Integer(1), sp.Integer(0)), (2, sp.Integer(1), sp.Integer(0))]),
        (sp.Integer(0),
         [(1, -sp.Integer(1) / n1, sp.Integer(1)),
          (2, -sp.Integer(1) / n2, sp.Integer(1))]),
        (-n3, [(3, sp.Integer(1), sp.Integer(0))]),
        (sp.Integer(0), [(3, -sp.Integer(1) / n3, sp.Integer(1))]),
    ]
    M, K = sp.zeros(4, 4), sp.zeros(4, 4)
    for i in range(4):
        for j in range(4):
            g, a = bilinear(gens[i], gens[j], sizes)
            K[i, j] = g
            M[i, j] = a
    return M, K


def numeric_bound(M, K, lval):
    Mn = np.array(sp.Matrix(M).subs(l, lval).evalf(), dtype=float)
    Kn = np.array(sp.Matrix(K).subs(l, lval).evalf(), dtype=float)
    ev, U = np.linalg.eigh(Kn)
    Ki = U @ np.diag(1.0 / np.sqrt(ev)) @ U.T
    S = Ki @ Mn @ Ki
    return float(np.linalg.eigvalsh((S + S.T) / 2).max())


def main():
    M, K = pencil(5, 5, l)
    print("M ="); print(M)
    print("\nK ="); print(K)

    spec = importlib.util.spec_from_file_location(
        "cb", "experiments/find_case_b_subspace.py")
    cb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cb)
    print("\ncross-validation against explicitly built graphs:")
    ok = True
    for lval in (5, 7, 11, 21, 101, 1001):
        N, A, offs, z, tips = cb.build((5, 5, lval))
        P, Q = cb.gens(N, offs, z, tips, (5, 5, lval))
        b1 = cb.top_of_span(A, [P[0] + P[1], Q[0] + Q[1], P[2], Q[2]])
        b2 = numeric_bound(M, K, lval)
        ok &= abs(b1 - b2) < 1e-9
        print(f"  l={lval:5d}: explicit {b1:.12f}  symbolic {b2:.12f}  "
              f"{'OK' if abs(b1-b2) < 1e-9 else 'MISMATCH'}")
    print("cross-validation:", "PASSED" if ok else "FAILED")

    D = sp.expand(M - TAU * K)
    print("\n2x2 principal minors of D = M - 2 sqrt2 K:")
    for i in range(4):
        for j in range(i + 1, 4):
            det = sp.factor(D[i, i] * D[j, j] - D[i, j] ** 2)
            print(f"  ({i},{j}): {det}")
    print("\ndiagonal of D:")
    for i in range(4):
        print(f"  D[{i},{i}] = {sp.factor(D[i, i])}")


if __name__ == "__main__":
    main()