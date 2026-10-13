"""Exploratory computation: signing norms of cubic graphs.

For a graph G and signing sigma: E(G) -> {+-1}, the signed adjacency matrix
A_sigma has (A_sigma)_{uv} = sigma({u,v}) if uv in E, else 0.

Define beta(G) = min_sigma ||A_sigma||  (operator norm = max |eigenvalue|).

Vertex switching sigma -> (delta_u sigma({u,v}) delta_v) preserves the spectrum
and leaves every sigma(v) unchanged. Switching along a spanning forest fixes the
signs of the forest edges, so representatives of switching classes are obtained
by taking the forest edges to +1 and freely choosing the remaining signs.
For a connected G this gives 2^(m-n+1) classes.
"""
from __future__ import annotations

import itertools
import numpy as np

RHO3 = 2 * np.sqrt(2)  # Ramanujan bound for d = 3


def adjacency(n, edges):
    A = np.zeros((n, n), dtype=np.int64)
    for u, v in edges:
        A[u, v] = 1
        A[v, u] = 1
    return A


def spanning_tree(n, edges):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    tree = []
    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
            tree.append((u, v))
    return tree


def signing_classes(n, edges):
    """Yield one representative signing per switching class (all as int arrays
    aligned with `edges`)."""
    tree = set(frozenset(e) for e in spanning_tree(n, edges))
    free_idx = [i for i, e in enumerate(edges) if frozenset(e) not in tree]
    mfree = len(free_idx)
    base = np.ones(len(edges), dtype=np.int64)
    for bits in itertools.product([1, -1], repeat=mfree):
        s = base.copy()
        s[free_idx] = bits
        yield s


def signed_matrix(n, edges, sigma):
    """sigma is indexed by edges (aligned with `edges`)."""
    M = np.zeros((n, n), dtype=np.int64)
    for i, (u, v) in enumerate(edges):
        M[u, v] = sigma[i]
        M[v, u] = sigma[i]
    return M.astype(float)


def beta_exact(n, edges, return_all=False):
    A = adjacency(n, edges)
    best, argbest = None, None
    hist = []
    for s in signing_classes(n, edges):
        As = signed_matrix(n, edges, s)
        w = np.linalg.eigvalsh(As)
        nrm = float(np.max(np.abs(w)))
        hist.append(nrm)
        if best is None or nrm < best - 1e-12:
            best, argbest = nrm, s.copy()
    return (best, argbest, hist) if return_all else (best, argbest)


def spectrum(n, edges):
    return np.sort(np.linalg.eigvalsh(adjacency(n, edges).astype(float)))[::-1]


def is_connected(n, edges):
    seen = {0}
    stack = [0]
    adj = {i: [] for i in range(n)}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == n


def is_cubic(n, edges):
    deg = [0] * n
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    return all(d == 3 for d in deg)


def ram_witness(n, edges):
    """Return (is_ramanujan, lambda_2) for a connected cubic graph."""
    w = spectrum(n, edges)
    return bool(w[1] <= RHO3 + 1e-9), float(w[1])


def has_perfect_matching(n, edges):
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    m = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)
    return 2 * len(m) == n


def cut_vertex16():
    edges = []
    base = 0
    centers = []
    for _ in range(3):
        a, b, c, d, z = base, base + 1, base + 2, base + 3, base + 4
        base += 5
        for (u, v) in [(a, c), (a, d), (b, c), (b, d), (c, d), (a, z), (b, z)]:
            edges.append((u, v))
        centers.append(z)
    for z in centers:
        edges.append((z, base))
    return base + 1, sorted(edges)


def k4():
    return 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def prism(n):
    edges = [(i, (i + 1) % n) for i in range(n)]
    edges += [(i, i + n) for i in range(n)]
    edges += [(i, (i + 1) % n + n) for i in range(n)]
    return 2 * n, sorted(edges)


def kneser_bipartite(a, b):
    """Cubic bipartite graph: cycle C on a+b alternating, plus a perfect matching."""
    n = a + b
    edges = [(i, (i + 1) % n) for i in range(n)]
    for i in range(a):
        edges.append((i, a + ((i * (b // a)) % b) if b % a == 0 else a + i))
    return n, sorted(set(edges))


CASES = [
    ("K4", k4()),
    ("triangular prism", prism(3)),
    ("cube", prism(4)),
    ("hexagonal prism", prism(6)),
]


def main():
    for name, (n, edges) in CASES:
        assert is_cubic(n, edges), name
        assert is_connected(n, edges), name
        ram, lam2 = ram_witness(n, edges)
        b, _ = beta_exact(n, edges)
        pm = has_perfect_matching(n, edges)
        print(f"{name:18s} n={n:3d} lam2={lam2:.6f} Ramanujan={ram} "
              f"beta={b:.6f} beta<={RHO3:.6f}:{b <= RHO3 + 1e-9} PM={pm}")

    print()
    n, edges = cut_vertex16()
    print("cut-vertex cubic n =", n, "cubic:", is_cubic(n, edges),
          "connected:", is_connected(n, edges))
    ram, lam2 = ram_witness(n, edges)
    print(f"  lam2={lam2:.6f} Ramanujan={ram} perfect matching={has_perfect_matching(n, edges)}")
    b, _ = beta_exact(n, edges)
    print(f"  beta={b:.6f}   (2*sqrt2 = {RHO3:.6f})  beta<=2sqrt2: {b <= RHO3 + 1e-9}")


if __name__ == "__main__":
    main()