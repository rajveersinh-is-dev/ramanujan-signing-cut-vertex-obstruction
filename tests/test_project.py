"""Tests for the signing/spectral computations and for the theorem steps.

Run:  python -m pytest tests -q      (or: python tests/test_project.py)
"""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from signings.core import (  # noqa: E402
    adjacency_matrix,
    beta,
    connected_components,
    degrees,
    girth,
    is_connected,
    is_cubic,
    is_ramanujan,
    min_lambda_max,
    min_lambda_min,
    perfect_matching,
    ram_bound,
    signed_adjacency_matrix,
    signing_classes,
    spectrum,
    vertex_sums,
)
from signings.theorems import (  # noqa: E402
    RHO3,
    WITNESS_L5,
    bottleneck,
    bottleneck_max_arms,
    det_02_exceeds_tau,
    gen_top,
    lambda2_of,
    psi,
    residual_pencil,
    subspace_bound,
    three_arm_graph,
    witness_gap,
)


# --------------------------------------------------------------------------- #
def k4():
    return 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def kn3():
    return 6, [(i, j) for i in range(3) for j in range(3, 6)]


def prism(k):
    E = [(i, (i + 1) % k) for i in range(k)]
    E += [(i + k, ((i + 1) % k) + k) for i in range(k)]
    E += [(i, i + k) for i in range(k)]
    return 2 * k, sorted(set(E))


def test_signing_classes_count():
    """Switching-class count is 2^(m-n+1) for connected G."""
    for (n, e) in (k4(), kn3(), prism(3), prism(5)):
        assert is_connected(n, e)
        nclasses = sum(1 for _ in signing_classes(n, e))
        assert nclasses == 2 ** (len(e) - n + 1), (n, nclasses)


def test_switching_invariance():
    """Vertex switching preserves the spectrum of A_sigma.

    (It does NOT preserve the vertex sums; that claim was made in an early
    draft of this file and is recorded as a corrected error in
    research/proof_attempts.md.)
    """
    rng = np.random.default_rng(0)
    n, e = prism(4)
    for _ in range(20):
        s = np.array([1 if rng.integers(2) else -1 for _ in e], dtype=np.int8)
        d = np.array([1 if rng.integers(2) else -1 for _ in range(n)], dtype=np.int64)
        sw = np.array([s[i] * d[u] * d[v] for i, (u, v) in enumerate(e)],
                      dtype=np.int8)
        w1 = np.sort(np.linalg.eigvalsh(signed_adjacency_matrix(n, e, s)))
        w2 = np.sort(np.linalg.eigvalsh(signed_adjacency_matrix(n, e, sw)))
        assert np.allclose(w1, w2, atol=1e-9)


def test_known_betas():
    n, e = k4()
    b, _ = beta(n, e)
    assert abs(b - np.sqrt(5)) < 1e-9, b
    n, e = kn3()
    b, _ = beta(n, e)
    assert abs(b - 2.0) < 1e-9, b
    n, e = prism(4)
    b, _ = beta(n, e)
    assert abs(b - np.sqrt(3)) < 1e-9, b


def test_beta_never_exceeds_d():
    for (n, e) in (k4(), kn3(), prism(3), prism(4), prism(5), prism(6)):
        b, _ = beta(n, e)
        assert b <= 3.0 + 1e-9


def test_mss_one_sided_holds_on_examples():
    """Marcus--Spielman--Srivastava: min_sigma lambda_max(A_sigma) <= 2 sqrt(d-1)."""
    for (n, e) in (k4(), kn3(), prism(3), prism(4), prism(5), prism(6)):
        assert min_lambda_max(n, e) <= RHO3 + 1e-9


def test_ram_bound():
    assert abs(ram_bound(3) - 2 * np.sqrt(2)) < 1e-12
    assert abs(ram_bound(4) - 2 * np.sqrt(3)) < 1e-12


def test_ramanujan_flags():
    n, e = prism(4)
    assert is_ramanujan(n, e, 3)
    n, e = three_arm_graph((5, 5, 5))
    assert not is_ramanujan(n, e, 3)


def test_girth():
    n, e = prism(4)
    assert girth(n, e) == 3
    n, e = three_arm_graph((5, 5, 5))
    assert girth(n, e) == 3


def test_perfect_matching():
    n, e = kn3()
    assert perfect_matching(n, e) is not None


# --------------------------------------------------------------------------- #
#  Theorem 1: the bottleneck bound
# --------------------------------------------------------------------------- #
def test_psi_bound():
    for r in np.linspace(1, 50, 400):
        assert psi(r) <= 1.0 + 1e-12
    assert abs(psi(1.0) - 1.0) < 1e-12


def test_bottleneck_matches_direct_computation():
    """3 - (a^2+b^2)/(ab(a+b)) is the Rayleigh quotient of
    (1/a) 1_{C_a} - (1/b) 1_{C_b} on the three-arm graph with arms a, b, a."""
    for (a, b) in [(5, 5), (5, 7), (7, 9), (5, 23), (11, 13)]:
        n, e = three_arm_graph((a, b, a))
        lam2 = lambda2_of(n, e)
        assert bottleneck(a, b) <= lam2 + 1e-9, (a, b)


def test_bottleneck_excludes_when_smaller_arm_ge_7():
    for (a, b) in [(7, 7), (7, 11), (9, 9), (11, 23), (99, 101)]:
        assert bottleneck(a, b) >= 20.0 / 7.0 - 1e-12
        assert bottleneck(a, b) > RHO3


# --------------------------------------------------------------------------- #
#  Theorem 2: a cubic graph with a cut vertex is never Ramanujan
# --------------------------------------------------------------------------- #
def test_residual_pencil_matches_explicit_vectors():
    for l in (5, 7, 9, 11, 21, 101):
        b1 = subspace_bound((5, 5, l), "residual")
        M, K = residual_pencil(l)
        b2 = gen_top(M, K)
        assert abs(b1 - b2) < 1e-9, (l, b1, b2)


def test_residual_bound_exceeds_rho3():
    for l in range(5, 200, 2):
        M, K = residual_pencil(l)
        assert gen_top(M, K) > RHO3, l


def test_det_02_negative_for_odd_l_ge_7():
    for l in range(7, 400, 2):
        assert det_02_exceeds_tau(l) < 0, l


def test_witness_gap_positive_at_l_5():
    gap = witness_gap(5)
    assert gap > 0, gap
    # exact value: 19584/625 - 2748 sqrt2 / 125
    expect = Fraction(19584, 625) - 2748 * np.sqrt(2) / 125
    assert abs(gap - float(expect)) < 1e-9


def test_three_arm_graphs_are_cubic_connected_with_cut_vertex():
    """G - z has exactly three components, one per prescribed arm order."""
    from signings.core import connected_components as cc
    for orders in [(5, 5, 5), (5, 7, 9), (5, 11, 13), (7, 7, 7), (9, 21, 23)]:
        n, e = three_arm_graph(orders)
        assert is_cubic(n, e)
        assert is_connected(n, e)
        z = n - 1
        remaining = [(u, v) for u, v in e if u != z and v != z]
        comps = cc(n - 1, [(u, v) for u, v in remaining if max(u, v) < z])
        assert len(comps) == 3, (orders, len(comps))
        assert sorted(len(c) for c in comps) == sorted(orders)
        # the centre has degree three and is the only point of attachment
        assert sum(1 for u, v in e if u == z or v == z) == 3


def test_no_three_arm_example_is_ramanujan():
    for orders in [(5, 5, 5), (5, 5, 7), (5, 5, 9), (5, 5, 11), (5, 7, 9),
                   (5, 9, 9), (7, 7, 7), (9, 11, 13), (5, 23, 5), (5, 5, 101)]:
        n, e = three_arm_graph(orders)
        lam2 = lambda2_of(n, e)
        assert lam2 > RHO3, (orders, lam2)


def _main() -> int:
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for f in fns:
        try:
            f()
            print(f"PASS  {f.__name__}")
        except Exception as exc:  # pragma: no cover
            failed += 1
            print(f"FAIL  {f.__name__}: {exc}")
    print(f"\n{len(fns) - failed}/{len(fns)} tests passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_main())