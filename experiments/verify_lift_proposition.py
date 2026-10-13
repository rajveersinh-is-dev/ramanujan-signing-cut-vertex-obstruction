"""Verification of the 2-lift proposition (Proposition 4 of proofs.md).

Builds explicit connected 2-lifts of cubic graphs G_0 from a non-constant
voltage g, checks:

  (i)  Spec(A(G)) = Spec(A(G_0)) sqcup Spec(A_{sigma_0}(G_0))   (the 2-lift
       spectral identity used in Proposition 4);
  (ii) if G is Ramanujan then ||A_{sigma_0}(G_0)|| <= 2 sqrt 2.

Run:  python experiments/verify_lift_proposition.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "src")
from signings.core import adjacency_matrix, is_cubic, is_connected, spectrum  # noqa: E402

RHO3 = 2 * np.sqrt(2)


def two_lift(n0, edges, s0):
    """The 2-lift of G_0 associated to the signing s0 (indexed by `edges`).

    Edge {u,v} with sign +1 lifts to {u^0 v^0, u^1 v^1}; with sign -1 it lifts
    to {u^0 v^1, u^1 v^0}.  Vertices are encoded as u + s*n0.
    """
    n = 2 * n0
    E = []
    for i, (u, v) in enumerate(edges):
        if s0[i] == 1:
            E += [(u, v), (u + n0, v + n0)]
        else:
            E += [(u, v + n0), (u + n0, v)]
    return n, sorted(set(E))


def signed_spectrum(n0, edges, s):
    A = np.zeros((n0, n0))
    for i, (u, v) in enumerate(edges):
        A[u, v] = s[i]
        A[v, u] = s[i]
    return np.linalg.eigvalsh(A)


def is_bipartite(n, edges):
    import collections
    adj = collections.defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    col = {}
    for s in range(n):
        if s in col:
            continue
        col[s] = 0
        st = [s]
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in col:
                    col[y] = 1 - col[x]
                    st.append(y)
                elif col[y] == col[x]:
                    return False
    return True


def main():
    bases = [
        ("K_4", 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]),
        ("K_{3,3}", 6, [(i, j) for i in range(3) for j in range(3, 6)]),
        ("prism C_3 x K_2", 6, [(0, 1), (1, 2), (2, 0),
                                (3, 4), (4, 5), (5, 3),
                                (0, 3), (1, 4), (2, 5)]),
        ("cube", 8, [(0, 1), (1, 2), (2, 3), (3, 0),
                     (4, 5), (5, 6), (6, 7), (7, 4),
                     (0, 4), (1, 5), (2, 6), (3, 7)]),
        ("Petersen", 10, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
                           (5, 7), (7, 9), (9, 6), (6, 8), (8, 5),
                           (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]),
    ]
    print(f"2*sqrt(2) = {RHO3:.12f}")
    print("For every signing s0 of G_0 we build the associated 2-lift G and check\n"
          "  (i)  Spec(A(G)) = Spec(A(G_0)) sqcup Spec(A_s0(G_0))  (2-lift identity);")
    print("  (ii) if G is Ramanujan then every eigenvalue of A_s0(G_0) lies in")
    print("       [-2sqrt2, 2sqrt2] except possibly the single value -3;")
    print("  (iii) if G_0 is BIPARTITE then no exceptional value occurs.")
    print("  (iv)  +3 never occurs in Spec(A_s0(G_0)).\n")
    ok = True
    tested = 0
    trivial = 0
    ram_lifts = 0
    ident_fail = 0
    exc_fail = 0
    bipartite_fail = 0
    pos3_fail = 0
    for (name, n0, e) in bases:
        assert is_cubic(n0, e), name
        g0_bipartite = is_bipartite(n0, e)
        m0 = len(e)
        A0 = np.sort(np.linalg.eigvalsh(adjacency_matrix(n0, e).astype(float)))
        shown = 0
        for mask in range(2 ** m0):
            s0 = np.array([1 if (mask >> i) & 1 else -1 for i in range(m0)],
                          dtype=np.int64)
            n, E = two_lift(n0, e, s0)
            assert is_cubic(n, E), name
            S0 = np.sort(signed_spectrum(n0, e, s0))
            AG = np.sort(spectrum(n, E))
            if not np.allclose(np.sort(np.concatenate([A0, S0])), AG,
                               atol=1e-8):
                ident_fail += 1
                ok = False
            if not is_connected(n, E):
                trivial += 1
                continue
            tested += 1
            rho = 2 * np.sqrt(2)
            nontriv = [x for x in AG if abs(abs(x) - 3.0) > 1e-9]
            ram = all(abs(x) <= rho + 1e-9 for x in nontriv)
            norm0 = float(np.max(np.abs(S0)))
            exceptional = norm0 > rho + 1e-9
            # (ii) at most one exceptional value, and it is -3
            if ram:
                ram_lifts += 1
                exc = [x for x in S0 if abs(x) > rho + 1e-9]
                if len(exc) > 1 or (exc and abs(exc[0] + 3.0) > 1e-8):
                    exc_fail += 1
                    ok = False
                # (iii) bipartite base => no exceptional value
                if g0_bipartite and exc:
                    bipartite_fail += 1
                    ok = False
            # (iv) +3 never occurs
            if any(abs(x - 3.0) < 1e-8 for x in S0):
                pos3_fail += 1
                ok = False
            if shown < 2:
                print(f"  {name:>18} n={n:3d}  G_0 bipartite={str(g0_bipartite):>5}"
                      f"  Ramanujan(G)={str(ram):>5}  ||A_s0(G_0)||={norm0:9.6f}"
                      f"  exceptional={bool(exc)}")
                shown += 1
    print()
    print(f"  signings of the bases examined: {tested + trivial} "
          f"(of which {trivial} give the trivial disconnected lift "
          f"G_0 x K_2)")
    print(f"  connected lifts: {tested}; of these Ramanujan: {ram_lifts}")
    print(f"  (i)   2-lift identity failures:               {ident_fail}")
    print(f"  (ii)  >1 exceptional value, or the value is not -3:  {exc_fail}")
    print(f"  (iii) bipartite base with an exceptional value present: {bipartite_fail}")
    print(f"  (iv)  +3 occurs in Spec(A_s0(G_0)):                  {pos3_fail}")
    print()
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
