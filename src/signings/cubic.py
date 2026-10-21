"""Exhaustive generation of connected cubic (3-regular) graphs.

Generation strategy.  Every cubic graph has a vertex; relabel so that this
vertex is 0 and its three neighbours are 1, 2, 3.  Then process the vertices
in the order 0, 1, ..., n-1; when vertex v is processed we must attach
(3 - deg(v)) further edges to vertices with larger index.  Vertices <= v are
"finalised" afterwards, so all still-open degree deficits must be realisable
inside the remaining vertex set -- this gives a cheap and *exact* feasibility
test.  Duplicates are removed with an isomorphism test (igraph canonical
permutation).

Verified against the classical counts of connected cubic graphs
(A002851): 1, 2, 5, 19, 85, 509, 4060 on 4, 6, 8, 10, 12, 14, 16 vertices.
"""
from __future__ import annotations

import itertools
from typing import Dict, List, Sequence, Set, Tuple

import numpy as np

from .core import Edge, is_connected

_igraph = None


def _ig():
    global _igraph
    if _igraph is None:
        import igraph
        _igraph = igraph
    return _igraph


def _canon(n: int, adj: Sequence[Set[int]]) -> bytes:
    """Canonical form of a graph: the adjacency matrix read off after applying
    igraph's canonical permutation.  Two graphs are isomorphic iff these agree."""
    ig = _ig()
    g = ig.Graph(n)
    for u in range(n):
        for v in adj[u]:
            if v > u:
                g.add_edge(u, v)
    perm = g.canonical_permutation()
    g_canon = g.permute_vertices(perm)
    return np.array(g_canon.get_adjacency().data, dtype=np.uint8).tobytes()


def _edges(adj: Sequence[Set[int]]) -> List[Edge]:
    e = []
    for u in range(len(adj)):
        for v in adj[u]:
            if v > u:
                e.append((u, v))
    return sorted(e)


def _feasible(adj: Sequence[Set[int]], n: int, v: int) -> bool:
    """After finalising vertex v, can the deficits on {v+1,...,n-1} be met?"""
    r = n - 1 - v  # number of not-yet-processed vertices
    if r == 0:
        return all(len(a) == 3 for a in adj)
    total = 0
    for u in range(v + 1, n):
        du = 3 - len(adj[u])
        if du < 0 or du > r - 1:
            return False
        total += du
    if total % 2 != 0:
        return False
    if total > r * (r - 1):
        return False
    return True


def generate_cubic(n: int, connected_only: bool = True) -> List[Tuple[int, List[Edge]]]:
    """All cubic graphs on n vertices, up to isomorphism, as (n, sorted edges)."""
    if n % 2 != 0 or n < 4:
        raise ValueError("cubic graphs need an even number of vertices >= 4")
    adj: List[Set[int]] = [set() for _ in range(n)]
    seen: Set[bytes] = set()
    out: List[Tuple[int, List[Edge]]] = []

    def rec(v: int) -> None:
        if v == n:
            if any(len(a) != 3 for a in adj):
                return
            if connected_only and not is_connected(n, _edges(adj)):
                return
            c = _canon(n, adj)
            if c in seen:
                return
            seen.add(c)
            out.append((n, _edges(adj)))
            return
        need = 3 - len(adj[v])
        rest = list(range(v + 1, n))
        if need < 0 or need > len(rest):
            return
        if need == 0:
            if _feasible(adj, n, v):
                rec(v + 1)
            return
        for S in itertools.combinations(rest, need):
            if any(x in adj[v] for x in S):
                continue
            ok = True
            for x in S:
                if len(adj[x]) >= 3:
                    ok = False
                    break
            if not ok:
                continue
            for x in S:
                adj[v].add(x)
                adj[x].add(v)
            if _feasible(adj, n, v):
                rec(v + 1)
            for x in S:
                adj[v].discard(x)
                adj[x].discard(v)

    for u in (1, 2, 3):
        adj[0].add(u)
        adj[u].add(0)
    rec(1)
    out.sort(key=lambda t: t[1])
    return out


def known_counts() -> Dict[int, int]:
    """Numbers of connected cubic graphs on n vertices (OEIS A002851)."""
    return {4: 1, 6: 2, 8: 5, 10: 19, 12: 85, 14: 509, 16: 4060,
            18: 41301, 20: 510489}
