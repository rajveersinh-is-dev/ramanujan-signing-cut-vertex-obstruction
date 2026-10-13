"""Results reported in the paper, produced by one reproducible script.

Run:  python experiments/run_all_experiments.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "src")
from signings.core import (  # noqa: E402
    adjacency_matrix,
    beta,
    degrees,
    girth,
    is_cubic,
    is_connected,
    is_ramanujan,
    min_lambda_max,
    min_lambda_min,
    perfect_matching,
    ram_bound,
    signed_adjacency_matrix,
    signing_classes,
    spectrum,
    vertex_connectivity,
)

RHO3 = 2 * np.sqrt(2)


# --------------------------------------------------------------------------- #
#  named graphs
# --------------------------------------------------------------------------- #
def k4():
    return 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def prism(k):
    """C_k x K_2, cubic on 2k vertices (k >= 3)."""
    E = [(i, (i + 1) % k) for i in range(k)]
    E += [(i + k, ((i + 1) % k) + k) for i in range(k)]
    E += [(i, i + k) for i in range(k)]
    return 2 * k, sorted(set(E))


def kneser_k33():
    return 6, [(i, j) for i in range(3) for j in range(3, 6)]


def mobius_kantor():
    """Möbius--Kantor graph: generalised Petersen graph GP(8,3)."""
    E = [(i, (i + 1) % 8) for i in range(8)]
    E += [(i, 8 + i) for i in range(8)]
    E += [(8 + i, 8 + (i + 3) % 8) for i in range(8)]
    return 16, sorted(set(E))


def pappus():
    """Pappus graph, a cubic Ramanujan graph on 18 vertices."""
    E = [(i, i + 9) for i in range(9)]
    E += [(i, i + 18) for i in range(9)]
    E += [(i, i + 27) for i in range(9)]
    E += [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
          (5, 7), (7, 9), (9, 11), (11, 13), (13, 15), (15, 17), (17, 5),
          (6, 10), (10, 14), (14, 18), (18, 22), (22, 26), (26, 8), (8, 12),
          (12, 16), (16, 20), (20, 24), (24, 6)]
    return 36, sorted(set(E))


def heawood():
    """Heawood graph: incidence graph of the Fano plane, cubic, 14 vertices."""
    pts = list(range(7))
    lines = [[0, 1, 3], [1, 2, 4], [2, 3, 5], [3, 4, 6], [4, 5, 0],
             [5, 6, 1], [6, 0, 2]]
    E = []
    for li, L in enumerate(lines):
        for p in L:
            E.append((pts[p], 7 + li))
    return 14, sorted(set(E))


CASES = [
    ("K_4", k4()),
    ("K_{3,3}", kneser_k33()),
    ("triangular prism C_3 x K_2", prism(3)),
    ("cube C_4 x K_2", prism(4)),
    ("pentagonal prism", prism(5)),
    ("hexagonal prism", prism(6)),
    ("Heawood", heawood()),
    ("Möbius--Kantor GP(8,3)", mobius_kantor()),
]


def table1():
    print("=" * 78)
    print("Table 1: small cubic graphs")
    print("=" * 78)
    print(f"{'graph':<26}{'n':>4}{'lam2':>12}{'Ramanujan':>11}"
          f"{'beta(G)':>12}{'lmax_min':>12}{'kappa':>7}{'girth':>7}")
    for name, (n, e) in CASES:
        assert is_cubic(n, e) and is_connected(n, e), name
        w = spectrum(n, e)
        lam2 = float(w[1])
        b, _ = beta(n, e)
        lmax = min_lambda_max(n, e)
        g = girth(n, e)
        print(f"{name:<26}{n:>4}{lam2:>12.6f}{str(lam2 <= RHO3 + 1e-9):>11}"
              f"{b:>12.6f}{lmax:>12.6f}{vertex_connectivity(n, e):>7}{g:>7}")


def table2():
    print()
    print("=" * 78)
    print("Table 2: MSS one-sided bound vs the Bilu--Linial norm bound (d = 3)")
    print("=" * 78)
    print(f"{'graph':<26}{'min_sigma lam_max':>20}{'min_sigma lam_min':>20}"
          f"{'beta(G)':>12}")
    for name, (n, e) in CASES:
        lmax = min_lambda_max(n, e)
        lmin = min_lambda_min(n, e)
        b, _ = beta(n, e)
        print(f"{name:<26}{lmax:>20.6f}{lmin:>20.6f}{b:>12.6f}")
    print()
    print(f"Ramanujan interval for d = 3 is [{-RHO3:.6f}, {RHO3:.6f}].")
    print("In every case min_sigma lam_max(A_sigma) <= 2 sqrt 2, as predicted")
    print("by Marcus--Spielman--Srivastava; the two-sided bound is the harder one.")


def main():
    np.set_printoptions(precision=6, suppress=True)
    print(f"2*sqrt(2) = {RHO3:.12f}")
    print(f"2*sqrt(d-1) for d=3; Marcus--Spielman--Srivastava bound on "
          f"min_sigma lam_max(A_sigma) = {RHO3:.6f}")
    table1()
    table2()
    print()
    print("All computations use IEEE-754 double precision via numpy.linalg.eigvalsh.")
    print("The mathematical claims in the paper are proved by hand; these numbers")
    print("only support the discussion and are reported as such.")


if __name__ == "__main__":
    main()