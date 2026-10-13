"""Exact witness for the residual case l = 5 of the cut-vertex theorem.

For odd l >= 7 the (P, R) 2x2 principal block of D = M - 2 sqrt(2) K has negative
determinant (hence a positive direction), which proves the bound exceeds
2 sqrt 2.  This script settles l = 5 by exhibiting a *rational* vector in the
4-dimensional span whose Rayleigh quotient exceeds 2 sqrt 2, and checking the
inequality exactly with sympy.
"""
from __future__ import annotations

import importlib.util
from fractions import Fraction

import numpy as np
import sympy as sp

l = sp.symbols('l', positive=True)
S2 = sp.sqrt(2)


def load():
    spec = importlib.util.spec_from_file_location(
        "cbs", "experiments/case_b_symbolic.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def find_witness(M, K, lval, grid=None):
    """Search for a small-denominator rational vector a with a'Ma > 2 sqrt2 a'Ka."""
    Mn = np.array(sp.Matrix(M).subs(l, lval).evalf(), dtype=float)
    Kn = np.array(sp.Matrix(K).subs(l, lval).evalf(), dtype=float)
    best = None
    if grid is None:
        grid = [Fraction(i, 5) for i in range(-10, 11)]
    import itertools
    for a in itertools.product(grid, repeat=4):
        if all(x == 0 for x in a):
            continue
        av = np.array([float(x) for x in a])
        num = av @ Mn @ av
        den = av @ Kn @ av
        if den <= 0:
            continue
        q = num / den
        if best is None or q > best[0]:
            best = (q, a)
    return best


def exact_check(M, K, lval, a):
    """Verify a'Ma - 2 sqrt2 a'Ka > 0 exactly, with rational a."""
    av = sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in a])
    Mv = sp.Matrix(M).subs(l, lval)
    Kv = sp.Matrix(K).subs(l, lval)
    expr = sp.expand((av.T * Mv * av)[0] - 2 * S2 * (av.T * Kv * av)[0])
    val = sp.N(expr, 30)
    return expr, sp.simplify(sp.radsimp(val) > 0), val


def main():
    m = load()
    M, K = m.pencil(5, 5, l)
    D = sp.expand(M - 2 * S2 * K)
    print("2x2 principal minor (P,R) of D = M - 2 sqrt2 K:")
    det02 = sp.factor(D[0, 0] * D[2, 2] - D[0, 2] ** 2)
    print("   ", det02)
    inner = sp.simplify(det02 / (-4 * (l + 11)))
    print("   factored:", sp.factor(inner))
    # inner > 0 iff l > (5 sqrt2 - 2)/(14 sqrt2 - 19)
    print("   threshold l >", sp.simplify((5 * S2 - 2) / (14 * S2 - 19)),
          "=", sp.N((5 * S2 - 2) / (14 * S2 - 19)))

    print("\ncase l = 5: searching for a rational witness")
    q, a = find_witness(M, K, 5)
    print(f"   best grid quotient = {q:.9f}  vector = {[str(x) for x in a]}")
    expr, ok, val = exact_check(M, K, 5, a)
    print("   exact expression a'Ma - 2 sqrt2 a'Ka =", expr)
    print("   numerical value:", val)
    print("   strictly positive:", ok)

    # also confirm the full symbolic bound at l=5 exceeds 2 sqrt 2 numerically
    b = m.numeric_bound(M, K, 5)
    print(f"   full 4-dim bound at l=5: {b:.12f}  > 2 sqrt2 = {2*np.sqrt(2):.12f}")


if __name__ == "__main__":
    main()