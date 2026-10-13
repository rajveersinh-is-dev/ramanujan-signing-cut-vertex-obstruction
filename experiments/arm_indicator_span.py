"""Optimal Rayleigh quotient on the span of the three arm indicators.

Let G be a connected cubic graph, z a cut vertex, and
G - z = C_1 U C_2 U C_3 with |C_i| = n_i (odd, >= 5); each C_i meets z in
exactly one edge, say at t_i.  Put S_i = C_i U {z} and A := a_1+a_2+a_3.

A vector in  W := span{ 1_{S_1}, 1_{S_2}, 1_{S_3} }  is
x = sum_i a_i 1_{S_i}, with x|_{C_i} = a_i and x_z = A.  Then

    ||x||^2 = sum_i a_i^2 n_i + A^2,
    x' A x  = sum_i a_i^2 (3 n_i - 1) + 2 A sum_i a_i
            = sum_i a_i^2 (3 n_i - 1) + 2 A^2,

and  <x,1> = 0  iff  sum_i a_i (n_i + 1) = 0.

So lambda_2(A(G)) >= R(n_1,n_2,n_3), the largest generalised eigenvalue of
(M, K) restricted to the 2-dimensional space c^perp, where

    M = diag(3 n_i - 1) + 2 1 1',   K = diag(n_i) + 1 1',   c = (n_i + 1).

Sanity check: the two-arm bottleneck vector (1/n_1) 1_{C_1} - (1/n_2) 1_{C_2}
lies in W, so R >= 3 - (n_1^2+n_2^2)/(n_1 n_2 (n_1+n_2)).
"""
from __future__ import annotations

import numpy as np

RHO3 = 2 * np.sqrt(2)


def _basis(c):
    c = c / np.linalg.norm(c)
    B = np.linalg.qr(np.column_stack([c, np.eye(3)]))[0]
    return B[:, 1:]


def best_quotient(n1, n2, n3):
    nn = np.array([n1, n2, n3], dtype=float)
    M = np.diag(3 * nn - 1) + 2 * np.ones((3, 3))
    K = np.diag(nn) + np.ones((3, 3))
    P = _basis(nn + 1.0)
    Mp, Kp = P.T @ M @ P, P.T @ K @ P
    ev, V = np.linalg.eigh(Kp)
    Ki = V @ np.diag(1.0 / np.sqrt(ev)) @ V.T   # Kp^{-1/2}
    S = Ki @ Mp @ Ki
    return float(np.linalg.eigvalsh((S + S.T) / 2).max())


def bottleneck(n1, n2):
    return 3.0 - (n1 * n1 + n2 * n2) / float(n1 * n2 * (n1 + n2))


def main():
    print(f"2*sqrt(2) = {RHO3:.9f}")
    b = bottleneck(5, 5)
    print(f"sanity: bottleneck(5,5) = {b:.9f}   R(5,5,5) = {best_quotient(5,5,5):.9f}"
          f"   (R must be >= bottleneck)")
    print()
    print(f"{'n1':>5} {'n2':>5} {'n3':>6} {'bottleneck':>12} {'R':>14}   >2sqrt2?")
    cases = [(5, 5, 5), (5, 5, 7), (5, 5, 9), (5, 5, 11), (5, 5, 101), (5, 5, 9999),
             (5, 7, 7), (5, 7, 9), (5, 7, 9999), (5, 9, 9), (7, 7, 7), (7, 9, 9999),
             (5, 11, 9999), (5, 21, 9999), (5, 51, 9999), (5, 101, 9999),
             (5, 9999, 9999), (9999, 9999, 9999), (5, 5, 101), (5, 5, 10001)]
    for (n1, n2, n3) in cases:
        R = best_quotient(n1, n2, n3)
        print(f"{n1:5d} {n2:5d} {n3:6d} {bottleneck(n1, n2):12.9f} {R:14.9f}"
              f"   {'YES' if R > RHO3 else 'NO'}")

    print("\nscan over odd 5 <= n1 <= n2 <= n3 <= 300:")
    odds = list(range(5, 301, 2))
    worst, at = 1e9, None
    bad = []
    for i, n1 in enumerate(odds):
        for j in range(i, len(odds)):
            n2 = odds[j]
            for k in range(j, len(odds)):
                n3 = odds[k]
                R = best_quotient(n1, n2, n3)
                if R < worst:
                    worst, at = R, (n1, n2, n3)
                if R <= RHO3 + 1e-12:
                    bad.append((n1, n2, n3, R))
    print(f"  min R = {worst:.9f} at {at}")
    print(f"  # scanned triples with R <= 2 sqrt 2: {len(bad)}")
    for x in bad[:20]:
        print("   ", x)


if __name__ == "__main__":
    main()