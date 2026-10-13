"""Independent randomised verification of Theorem 2.

Theorem 2 says: no connected cubic graph with a cut vertex is Ramanujan, i.e.
lambda_2(A(G)) > 2 sqrt(2).

This script checks it on a large random sample of three-arm graphs built by
randomised gluing, and compares the *optimal* Rayleigh quotient found on the
full space {x : <x,1> = 0} (i.e. the true lambda_2) against 2 sqrt 2.  It also
checks that the *provable* lower bounds used in the proof (Corollary 1 for
n_2 >= 7, and the Theorem 3 pencil for n_1 = n_2 = 5) never exceed the true
lambda_2.

Run:  python experiments/verify_theorem2_random.py
"""
from __future__ import annotations

import random
import sys

import numpy as np

sys.path.insert(0, "src")
from signings.core import is_cubic, is_connected  # noqa: E402
from signings.theorems import (  # noqa: E402
    RHO3,
    bottleneck,
    gen_top,
    lambda2_of,
    residual_pencil,
    subspace_bound,
    three_arm_graph,
)

SEED = 20261013
SAMPLES = 4000


def random_odd(nmin: int, nmax: int, rng: random.Random) -> int:
    vals = list(range(nmin, nmax + 1, 2))
    return rng.choice(vals)


def random_arm(n: int, rng: random.Random):
    """A random arm of odd order n >= 5: a random cubic graph on n-1 vertices
    with one edge subdivided into a path of length 2 (so the new vertex has
    degree 2).  Falls back to the canonical arm when the subdivision would
    disconnect."""
    import itertools
    # random cubic graph on n-1 vertices via a simple randomised pairing
    m = n - 1
    if m < 4:
        from signings.theorems import arm_of_order
        return arm_of_order(n)
    for _ in range(60):
        # random 3-regular multigraph then simple-ise by rejection
        stubs = [(v, i) for v in range(m) for i in range(3)]
        rng.shuffle(stubs)
        edges = []
        ok = True
        used = set()
        for k in range(0, len(stubs), 2):
            (u, _), (v, _) = stubs[k], stubs[k + 1]
            if u == v or frozenset((u, v)) in used:
                ok = False
                break
            used.add(frozenset((u, v)))
            edges.append((u, v))
        if not ok or len(edges) != 3 * m // 2:
            continue
        # connected?
        import collections
        adj = collections.defaultdict(list)
        for u, v in edges:
            adj[u].append(v)
            adj[v].append(u)
        seen = {0}
        st = [0]
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    st.append(y)
        if len(seen) != m:
            continue
        # subdivide one edge
        (x, y) = edges[0]
        w = m
        E = [e for e in edges if frozenset(e) != frozenset((x, y))]
        E += [(x, w), (w, y)]
        return n, sorted(E), w
    from signings.theorems import arm_of_order
    return arm_of_order(n)


def random_three_arm_graph(nmax: int, rng: random.Random):
    off = 0
    edges, tips = [], []
    orders = []
    for _ in range(3):
        o = random_odd(5, nmax, rng)
        n, e, t = random_arm(o, rng)
        edges += [(u + off, v + off) for u, v in e]
        tips.append(off + t)
        orders.append(n)
        off += n
    z = off
    edges += [(t, z) for t in tips]
    return sorted(orders), z + 1, sorted(set(edges))


def main():
    rng = random.Random(SEED)
    nprng = np.random.default_rng(SEED)
    print(f"seed = {SEED}, samples = {SAMPLES}")
    print(f"2*sqrt(2) = {RHO3:.12f}\n")
    worst_margin = 1e9
    worst_at = None
    bad = []
    bound_violations = []
    for it in range(SAMPLES):
        orders, n, e = random_three_arm_graph(11, rng)
        orders.sort()
        assert is_cubic(n, e) and is_connected(n, e)
        lam2 = lambda2_of(n, e)
        margin = lam2 - RHO3
        if margin < worst_margin:
            worst_margin, worst_at = margin, (orders, n, lam2)
        if lam2 <= RHO3 + 1e-9:
            bad.append((orders, n, lam2))
        # the provable bound used in the proof must never exceed lambda_2
        if orders[1] >= 7:
            bnd = max(bottleneck(orders[1], orders[2]),
                      bottleneck(orders[1], orders[2]))
        else:
            M, K = residual_pencil(orders[2])
            bnd = gen_top(M, K)
        if bnd > lam2 + 1e-8:
            bound_violations.append((orders, bnd, lam2))
    print(f"checked {SAMPLES} random three-arm cubic graphs (arm orders 5..11)")
    print(f"  graphs that ARE Ramanujan: {len(bad)}   (expected 0)")
    for b in bad[:10]:
        print("    ", b)
    print(f"  smallest margin lambda_2 - 2sqrt(2) = {worst_margin:.9f}")
    print(f"    attained at arm orders {worst_at[0]}, n = {worst_at[1]}, "
          f"lambda_2 = {worst_at[2]:.9f}")
    print(f"  proof-bound violations (bound > true lambda_2): "
          f"{len(bound_violations)}")
    for v in bound_violations[:10]:
        print("    ", v)
    ok = (len(bad) == 0) and (len(bound_violations) == 0)
    print()
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())