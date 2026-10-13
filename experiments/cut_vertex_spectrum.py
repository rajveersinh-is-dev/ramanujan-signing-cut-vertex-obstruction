"""Is every connected cubic graph with a cut vertex non-Ramanujan?

For a connected cubic graph G and a cut vertex z, G - z has exactly three
components C_1, C_2, C_3, each meeting z in exactly one edge, and each of odd
order >= 5.  Taking x = (1/|C_1|) 1_{C_1} - (1/|C_2|) 1_{C_2} (perpendicular to
the all-ones vector) gives

    lambda_2(A(G))  >=  3 - (n_1^2 + n_2^2) / (n_1 n_2 (n_1 + n_2)).

This script builds cubic graphs with a cut vertex having prescribed small arms
and computes the exact lambda_2, testing whether the bound excludes Ramanujan.
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "src")
from signings.core import is_cubic, is_connected, spectrum  # noqa: E402

RHO3 = 2 * np.sqrt(2)


def cubic_h(m):
    """Some cubic graph on m (even, m >= 4) vertices."""
    if m == 4:
        return 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    k = m // 2
    assert k >= 3
    E = [(i, (i + 1) % k) for i in range(k)]
    E += [(i + k, ((i + 1) % k) + k) for i in range(k)]
    E += [(i, i + k) for i in range(k)]
    return 2 * k, sorted(set(E))


def arm(order):
    """Component of odd order `order`, cubic except one vertex of degree 2."""
    assert order % 2 == 1 and order >= 5, order
    hn, he = cubic_h(order - 1)
    E = [e for e in he if frozenset(e) != frozenset((0, 1))]
    w = hn
    E = E + [(w, 0), (w, 1)]
    return order, sorted(E)


def degree2_vertex(n, edges):
    deg = [0] * n
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    w = [i for i in range(n) if deg[i] == 2]
    assert len(w) == 1, deg
    return w[0]


def three_arms(orders):
    """Cubic graph with cut vertex z whose three arms have the given orders."""
    off = 0
    edges = []
    tips = []
    for o in orders:
        n, e = arm(o)
        edges += [(u + off, v + off) for u, v in e]
        tips.append(off + degree2_vertex(n, e))
        off += n
    z = off
    edges += [(w, z) for w in tips]
    return z + 1, sorted(set(edges))


def bottleneck_bound(a, b):
    return 3.0 - (a * a + b * b) / float(a * b * (a + b))


def main():
    print(f"Ramanujan threshold for d = 3: 2*sqrt(2) = {RHO3:.9f}\n")
    print("bottleneck bound B(a,b) = 3 - (a^2+b^2)/(ab(a+b)), for a <= b odd, >= 5:")
    for (a, b) in [(5, 5), (5, 7), (5, 9), (5, 11), (5, 13), (5, 15),
                   (7, 7), (7, 9), (7, 13), (9, 9), (9, 13), (11, 11)]:
        B = bottleneck_bound(a, b)
        print(f"  a={a:2d} b={b:2d}  B={B:.9f}   "
              f"{'EXCLUDES Ramanujan' if B > RHO3 else 'does NOT exclude'}")

    print("\nexact lambda_2 for three-arm graphs (arms = a, b, c):")
    print(f"{'a':>3}{'b':>3}{'c':>4}{'n':>6}{'B(a,b)':>14}{'lam2':>16}  Ramanujan?")
    rows = []
    for (a, b) in [(5, 5), (5, 7), (5, 9), (5, 11), (5, 13), (5, 15), (5, 21),
                   (5, 31), (5, 51), (7, 7), (7, 9), (7, 13), (7, 21),
                   (9, 9), (9, 13), (9, 21), (11, 11), (11, 21)]:
        for c in (5, 7, 9, 11, 13, 15, 21, 31, 51, 101):
            if c < b:
                continue
            n, e = three_arms((a, b, c))
            assert is_cubic(n, e), (n, e)
            assert is_connected(n, e)
            w = spectrum(n, e)
            lam2 = float(w[1])
            rows.append((a, b, c, n, lam2, lam2 <= RHO3 + 1e-9))
            print(f"{a:3d}{b:3d}{c:4d}{n:6d}{bottleneck_bound(a,b):14.9f}"
                  f"{lam2:16.9f}  {lam2 <= RHO3 + 1e-9}")
    bad = [r for r in rows if r[5]]
    print()
    if bad:
        print("!! cubic graphs with a cut vertex that ARE Ramanujan:")
        for r in bad:
            print("   arms", r[:3], "n =", r[3], "lam2 =", r[4])
    else:
        print("Result: none of these three-arm graphs is Ramanujan.")


if __name__ == "__main__":
    main()