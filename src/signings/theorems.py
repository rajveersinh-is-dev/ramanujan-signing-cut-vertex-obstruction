"""Executable versions of the theorems proved in mathematics/proofs.md.

Every routine here is a direct transcription of a step of a hand proof.  They
are used by the test-suite and by experiments/run_all_experiments.py; the
mathematical content is the proofs in mathematics/proofs.md, not this code.
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from typing import Dict, List, Sequence, Tuple

import numpy as np

Edge = Tuple[int, int]

RHO3 = 2 * np.sqrt(2)


# --------------------------------------------------------------------------- #
# Theorem 1: the two-arm bottleneck bound.
# --------------------------------------------------------------------------- #
def psi(r: float) -> float:
    """psi(r) = (r + 1/r)/(1 + r); equals 1 at r = 1 and is <= 1 for r >= 1."""
    return (r + 1.0 / r) / (1.0 + r)


def bottleneck(a: int, b: int) -> float:
    """3 - (a^2 + b^2)/(ab(a+b)).  The bound of Theorem 1 for arms of
    orders a <= b."""
    return 3.0 - (a * a + b * b) / float(a * b * (a + b))


def bottleneck_max_arms(n: Sequence[int]) -> float:
    return max(bottleneck(n[i], n[j])
               for i, j in itertools.combinations(range(len(n)), 2))


def psi_bound(n: int, r: float) -> float:
    """3 - psi(r)/n, the form used in the proof (n = smaller arm order)."""
    return 3.0 - psi(r) / n


# --------------------------------------------------------------------------- #
# The four-dimensional test subspace for the residual case n_1 = n_2 = 5.
# --------------------------------------------------------------------------- #
def residual_pencil(l: int) -> Tuple[np.ndarray, np.ndarray]:
    """(M, K) for the span of

        P = p_1 + p_2,  Q = q_1 + q_2,  R = p_3,  S = q_3

    where p_i = 1_{C_i} - n_i e_z, q_i = e_{t_i} - (1/n_i) 1_{C_i}, n_1=n_2=5,
    n_3 = l.  Computed with exact rational arithmetic where possible."""
    L = Fraction(l)
    five = Fraction(5)
    M = np.array([
        [Fraction(-12), Fraction(-88, 5), -2 * L - 10, -10 + Fraction(10, 1) / L],
        [Fraction(-88, 5), Fraction(-12, 25), -Fraction(8, 5) * L, Fraction(0)],
        [-2 * L - 10, -Fraction(8, 5) * L, L - 1, -L + Fraction(1, 1) / L],
        [-10 + Fraction(10, 1) / L, Fraction(0), -L + Fraction(1, 1) / L,
         -Fraction(1, 1) / L - Fraction(1, 1) / (L * L)],
    ], dtype=object)
    K = np.array([
        [Fraction(110), Fraction(0), 10 * L, Fraction(0)],
        [Fraction(0), Fraction(8, 5), Fraction(0), Fraction(0)],
        [10 * L, Fraction(0), L * L + L, Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(0), 1 - Fraction(1, 1) / L],
    ], dtype=object)
    return M, K


def gen_top(M: np.ndarray, K: np.ndarray) -> float:
    """Largest generalised eigenvalue of the pencil (M, K), numerically."""
    Mn = np.array([[float(x) for x in row] for row in M], dtype=float)
    Kn = np.array([[float(x) for x in row] for row in K], dtype=float)
    ev, U = np.linalg.eigh(Kn)
    Ki = U @ np.diag(1.0 / np.sqrt(ev)) @ U.T
    S = Ki @ Mn @ Ki
    return float(np.linalg.eigvalsh((S + S.T) / 2).max())


def det_02_exceeds_tau(l: int) -> float:
    """Value of the (0,2) principal minor of M - 2 sqrt2 K, which is

        -4 (l + 11) ( (14 sqrt2 - 19) l - 5 sqrt2 + 2 ),

    and is < 0 for every odd l >= 7."""
    v = -4.0 * (l + 11) * ((14 * np.sqrt(2) - 19) * l - 5 * np.sqrt(2) + 2)
    return v


WITNESS_L5 = (Fraction(-3, 5), Fraction(1, 5), Fraction(6, 5), Fraction(-2, 5))


def witness_gap(l: int) -> float:
    """a'Ma - 2 sqrt2 a'Ka for the rational witness a = WITNESS_L5 at l = 5."""
    M, K = residual_pencil(l)
    a = np.array(WITNESS_L5, dtype=object)
    num = float(a @ M @ a)
    den = float(a @ K @ a)
    return num - 2 * np.sqrt(2) * den


# --------------------------------------------------------------------------- #
# Explicit construction of the three-arm graph (used for cross-checking).
# --------------------------------------------------------------------------- #
def arm_of_order(n: int) -> Tuple[int, List[Edge], int]:
    """An arm of odd order n: cubic except for the tip, which has degree 2."""
    assert n % 2 == 1 and n >= 5, n
    if n == 5:
        he = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    else:
        k = (n - 1) // 2
        he = sorted(set([(i, (i + 1) % k) for i in range(k)]
                        + [(i + k, ((i + 1) % k) + k) for i in range(k)]
                        + [(i, i + k) for i in range(k)]))
    E = [e for e in he if frozenset(e) != frozenset((0, 1))]
    return n, sorted(E + [(n - 1, 0), (n - 1, 1)]), n - 1


def three_arm_graph(orders: Sequence[int]) -> Tuple[int, List[Edge]]:
    """A connected cubic graph with a cut vertex z and arms of the given
    (odd, >= 5) orders."""
    off = 0
    edges: List[Edge] = []
    tips: List[int] = []
    for o in orders:
        nn, e, t = arm_of_order(o)
        edges += [(u + off, v + off) for u, v in e]
        tips.append(off + t)
        off += nn
    z = off
    edges += [(t, z) for t in tips]
    return z + 1, sorted(edges)


def lambda2_of(n: int, edges: Sequence[Edge]) -> float:
    A = np.zeros((n, n))
    for u, v in edges:
        A[u, v] = A[v, u] = 1.0
    return float(np.linalg.eigvalsh(A)[::-1][1])


def subspace_bound(n_orders: Sequence[int],
                   generators: str = "residual") -> float:
    """Top Rayleigh quotient on the test subspace, built from explicit
    coordinate vectors on the concrete graph `three_arm_graph(n_orders)`."""
    n, edges = three_arm_graph(n_orders)
    A = np.zeros((n, n))
    for u, v in edges:
        A[u, v] = A[v, u] = 1.0
    off, tips, offs = 0, [], []
    for o in n_orders:
        nn, _, t = arm_of_order(o)
        tips.append(off + t)
        offs.append(off)
        off += nn
    z = off
    P, Q = [], []
    for i, o in enumerate(n_orders):
        p = np.zeros(n)
        p[offs[i]:offs[i] + o] = 1
        p[z] = -o
        P.append(p)
        q = np.zeros(n)
        q[tips[i]] += 1
        q[offs[i]:offs[i] + o] -= 1.0 / o
        Q.append(q)
    if generators == "residual":
        V = [P[0] + P[1], Q[0] + Q[1], P[2], Q[2]]
    elif generators == "full6":
        V = P + Q
    elif generators == "bottleneck":
        V = [P[1], P[2]]
    else:
        raise ValueError(generators)
    M = np.column_stack(V)
    assert abs(M.sum(axis=0)).max() < 1e-8
    K = M.T @ M
    Aa = M.T @ (A @ M)
    ev, U = np.linalg.eigh(K)
    Ki = U @ np.diag(1.0 / np.sqrt(ev)) @ U.T
    S = Ki @ Aa @ Ki
    return float(np.linalg.eigvalsh((S + S.T) / 2).max())