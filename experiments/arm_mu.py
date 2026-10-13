"""Refined obstruction: mean-zero test vectors inside the arms of a cut vertex.

If G is a connected cubic graph and z is a cut vertex, then G - z has exactly
three components C_1, C_2, C_3, each joined to z by a single edge and each of
odd order >= 5.  Any vector x with x_z = 0 that is mean-zero on each C_i and
supported on V(G) \ {z} satisfies <x, 1> = 0, so

    x' A(G) x = sum_i (x|_{C_i})' A_{C_i} (x|_{C_i}).

Hence  lambda_2(A(G)) >= max_i mu(C_i),  where
mu(H) = max { u' A_H u / ||u||^2 : u != 0, <u, 1_H> = 0 }
is the largest mean-zero Rayleigh quotient of H.
"""
from __future__ import annotations

import itertools
import sys

import numpy as np

sys.path.insert(0, "src")
from signings.core import adjacency_matrix, spectrum  # noqa: E402

RHO3 = 2 * np.sqrt(2)


def mu(n, edges):
    """max of u'A u /||u||^2 over u != 0 with sum u = 0."""
    A = adjacency_matrix(n, edges).astype(float)
    n = A.shape[0]
    # projector onto mean-zero subspace
    P = np.eye(n) - np.ones((n, n)) / n
    B = P @ A @ P
    w = np.linalg.eigvalsh(B)
    return float(max(w.max(), -w.min())), A


def all_arms_upto(order):
    """All (n, edges) 'arm' graphs of odd order n: cubic except one degree-2
    vertex.  Built as: take any cubic graph H on n-1 vertices and an edge
    {x,y}; delete that edge and add a vertex w adjacent to x and y."""
    from signings.cubic import generate_cubic
    out = []
    for (hn, he) in generate_cubic(order - 1, connected_only=True):
        for (x, y) in he:
            E = [e for e in he if frozenset(e) != frozenset((x, y))]
            E = E + [(hn, x), (hn, y)]
            out.append((order, sorted(E)))
    # dedupe arms up to isomorphism of the whole arm
    import igraph as ig
    seen = set()
    uniq = []
    for (n, e) in out:
        g = ig.Graph(n)
        for u, v in e:
            g.add_edge(u, v)
        M = np.zeros((n, n), dtype=np.uint8)
        p = g.canonical_permutation()
        for u, v in e:
            M[p[u], p[v]] = 1
            M[p[v], p[u]] = 1
        b = M.tobytes()
        if b not in seen:
            seen.add(b)
            uniq.append((n, e))
    return uniq


def main():
    print(f"2*sqrt(2) = {RHO3:.9f}\n")
    for n in (5, 7, 9):
        arms = all_arms_upto(n)
        vals = []
        for (nn, e) in arms:
            m, _ = mu(nn, e)
            vals.append(m)
        vals = np.array(vals)
        print(f"order {n}: {len(arms)} arm types; "
              f"min mu = {vals.min():.9f}, max mu = {vals.max():.9f}, "
              f"{'ALL > 2sqrt2' if vals.min() > RHO3 else 'some <= 2sqrt2'}")


if __name__ == "__main__":
    main()
