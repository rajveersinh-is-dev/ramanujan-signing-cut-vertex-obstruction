"""Core objects for the study of signings of regular graphs.

A signing of a graph G is a map sigma: E(G) -> {-1, +1}.  The *signed adjacency
matrix* A_sigma(G) is the symmetric n x n matrix with

    (A_sigma)_{uv} = sigma({u,v})  if {u,v} in E(G),   0 otherwise,

with zero diagonal.  Its spectrum is invariant under *vertex switching*,
sigma -> (delta_u sigma({u,v}) delta_v)_{u,v} for any delta in {-1,+1}^V:
indeed A_switched = D A_sigma D with D = diag(delta).

The vertex sums sigma(v) := sum_{e ~ v} sigma(e) are NOT invariant under
switching (only their values after switching are sigma'(v) = delta_v
* sum_{w ~ v} sigma(uv) delta_w).  Two useful consequences of switching are
used below: the orbit of a signing has size 2^(n-1), and switching along a
spanning forest fixes the signs of the forest edges, which is how
`signing_classes` enumerates one representative per orbit.

Note also that max_v |sigma(v)| does *not* bound ||A_sigma||: for K_4 signed
with a perfect matching all vertex sums are +-1 while the spectral norm is 3.

We write  beta(G) := min_sigma ||A_sigma(G)||_op
with ||.||_op the Euclidean operator norm (max modulus of an eigenvalue, since
A_sigma is symmetric).
"""
from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

import numpy as np

Edge = Tuple[int, int]


def adjacency_matrix(n: int, edges: Sequence[Edge]) -> np.ndarray:
    A = np.zeros((n, n), dtype=np.int64)
    for u, v in edges:
        if u == v:
            raise ValueError("loops are not allowed")
        A[u, v] += 1
        A[v, u] += 1
    return A


def signed_adjacency_matrix(n: int, edges: Sequence[Edge],
                            sigma: Sequence[int]) -> np.ndarray:
    """sigma is aligned with `edges`."""
    if len(sigma) != len(edges):
        raise ValueError("signing length mismatch")
    M = np.zeros((n, n), dtype=np.int64)
    for i, (u, v) in enumerate(edges):
        if sigma[i] not in (-1, 1):
            raise ValueError("signing values must be +-1")
        M[u, v] = sigma[i]
        M[v, u] = sigma[i]
    return M


def spectrum(n: int, edges: Sequence[Edge]) -> np.ndarray:
    """Eigenvalues of the unsigned adjacency matrix, in decreasing order."""
    A = adjacency_matrix(n, edges).astype(float)
    return np.linalg.eigvalsh(A)[::-1].copy()


def vertex_sums(n: int, edges: Sequence[Edge], sigma: Sequence[int]) -> np.ndarray:
    s = np.zeros(n, dtype=np.int64)
    for i, (u, v) in enumerate(edges):
        s[u] += sigma[i]
        s[v] += sigma[i]
    return s


def degrees(n: int, edges: Sequence[Edge]) -> List[int]:
    deg = [0] * n
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    return deg


def is_regular(n: int, edges: Sequence[Edge], d: int | None = None) -> bool:
    deg = degrees(n, edges)
    return len(set(deg)) == 1 and (d is None or deg[0] == d)


def is_cubic(n: int, edges: Sequence[Edge]) -> bool:
    return is_regular(n, edges, 3)


def neighbour_sets(n: int, edges: Sequence[Edge]) -> List[set]:
    N: List[set] = [set() for _ in range(n)]
    for u, v in edges:
        N[u].add(v)
        N[v].add(u)
    return N


def connected_components(n: int, edges: Sequence[Edge]) -> List[List[int]]:
    N = neighbour_sets(n, edges)
    seen = [False] * n
    comps: List[List[int]] = []
    for s in range(n):
        if seen[s]:
            continue
        stack, comp = [s], []
        seen[s] = True
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in N[x]:
                if not seen[y]:
                    seen[y] = True
                    stack.append(y)
        comps.append(sorted(comp))
    return comps


def is_connected(n: int, edges: Sequence[Edge]) -> bool:
    return len(connected_components(n, edges)) == 1


def girth(n: int, edges: Sequence[Edge]) -> int:
    """Length of the shortest cycle; infinity if acyclic."""
    import heapq
    N = neighbour_sets(n, edges)
    INF = float("inf")
    best = INF
    for s in range(n):
        dist = {s: 0}
        par = {s: -1}
        pq = [(0, s)]
        while pq:
            d, x = heapq.heappop(pq)
            if d > dist.get(x, INF):
                continue
            if d >= best:
                break
            for y in N[x]:
                if y == par[x]:
                    continue
                nd = d + 1
                if y not in dist:
                    dist[y] = nd
                    par[y] = x
                    heapq.heappush(pq, (nd, y))
                else:
                    best = min(best, d + nd)
    return int(best) if best != INF else -1


def vertex_connectivity(n: int, edges: Sequence[Edge]) -> int:
    """Vertex connectivity k(G).  Uses networkx when available (exact); the
    brute-force fallback is only intended for very small graphs."""
    if n <= 1:
        return 0
    try:
        import networkx as nx
        G = nx.Graph()
        G.add_nodes_from(range(n))
        G.add_edges_from(edges)
        return int(nx.node_connectivity(G))
    except Exception:  # pragma: no cover - fallback
        return min(degrees(n, edges))


def _min_separating_set(n: int, edges: Sequence[Edge], u: int, v: int) -> int:
    """Smallest S (disjoint from {u,v}) whose removal separates u from v."""
    N = neighbour_sets(n, edges)
    others = [x for x in range(n) if x not in (u, v)]
    import itertools
    for r in range(0, n - 1):
        for S in itertools.combinations(others, r):
            blocked = set(S)
            seen = {u}
            stack = [u]
            while stack:
                x = stack.pop()
                for y in N[x]:
                    if y not in blocked and y not in seen:
                        seen.add(y)
                        stack.append(y)
            if v not in seen:
                return r
    return n


def perfect_matching(n: int, edges: Sequence[Edge]) -> List[Edge] | None:
    """Return a perfect matching, or None.  Blossom algorithm via networkx."""
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    M = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)
    M = sorted(tuple(sorted(e)) for e in M)
    if 2 * len(M) == n:
        return M
    return None


def signing_classes(n: int, edges: Sequence[Edge]) -> Iterable[np.ndarray]:
    """One representative per vertex-switching class.

    Switching along a spanning forest fixes the signs of the forest edges and
    leaves the remaining signs free; two signings in the same switching class
    are obtained from one another this way.  For a connected G on n vertices
    with m edges this yields exactly 2^(m-n+1) representatives.
    """
    tree = set(frozenset(e) for e in _spanning_forest(n, edges))
    free = [i for i, e in enumerate(edges) if frozenset(e) not in tree]
    mfree = len(free)
    import itertools
    for bits in itertools.product((1, -1), repeat=mfree):
        s = np.ones(len(edges), dtype=np.int8)
        if free:
            s[free] = bits
        yield s


def _spanning_forest(n: int, edges: Sequence[Edge]) -> List[Edge]:
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    out: List[Edge] = []
    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
            out.append((u, v))
    return out


def beta(n: int, edges: Sequence[Edge],
         return_all: bool = False):
    """Exact beta(G) by exhaustive search over switching classes."""
    norms = []
    best, argbest = None, None
    for s in signing_classes(n, edges):
        M = signed_adjacency_matrix(n, edges, s).astype(float)
        w = np.linalg.eigvalsh(M)
        nrm = float(np.max(np.abs(w)))
        norms.append(nrm)
        if best is None or nrm < best - 1e-12:
            best, argbest = nrm, s.copy()
    if return_all:
        return best, argbest, norms
    return best, argbest


def min_lambda_max(n: int, edges: Sequence[Edge]) -> float:
    """min_sigma lambda_max(A_sigma) -- the quantity bounded by MSS."""
    best = None
    for s in signing_classes(n, edges):
        M = signed_adjacency_matrix(n, edges, s).astype(float)
        v = float(np.max(np.linalg.eigvalsh(M)))
        if best is None or v < best - 1e-12:
            best = v
    return best


def min_lambda_min(n: int, edges: Sequence[Edge]) -> float:
    """max_sigma lambda_min(A_sigma) (equivalently -min_sigma lambda_max(A_{-sigma}))."""
    best = None
    for s in signing_classes(n, edges):
        M = signed_adjacency_matrix(n, edges, s).astype(float)
        v = float(np.min(np.linalg.eigvalsh(M)))
        if best is None or v > best + 1e-12:
            best = v
    return best


def ram_bound(d: int) -> float:
    return 2.0 * np.sqrt(d - 1)


def is_ramanujan(n: int, edges: Sequence[Edge], d: int) -> bool:
    """Connected d-regular G is Ramanujan iff every eigenvalue except the
    trivial ones (+/- d) lies in [-2 sqrt(d-1), 2 sqrt(d-1)]."""
    w = spectrum(n, edges)
    rho = ram_bound(d)
    # trivial: d, and -d if bipartite
    triv = [d]
    if w[-1] < -d + 1e-9:
        triv.append(-d)
    nontriv = [x for x in w if abs(x - d) > 1e-9 and abs(x + d) > 1e-9]
    return all(abs(x) <= rho + 1e-9 for x in nontriv)