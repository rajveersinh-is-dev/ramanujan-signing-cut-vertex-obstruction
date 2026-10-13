'''Optimal Rayleigh quotient on the span of the three arm indicators.

Let G be a connected cubic graph, z a cut vertex, and
G - z = C_1 \cup C_2 \cup C_3 with |C_i| = n_i (odd, >= 5); each C_i meets z in
exactly one edge, say at t_i.  Put S_i = C_i \cup {z} and A := a_1+a_2+a_3.

A vector in  W := span{ 1_{S_1}, 1_{S_2}, 1_{S_3} }  is
x = sum_i a_i 1_{S_i}, with x|_{C_i} = a_i and x_z = A.  Then

    ||x||^2 = sum_i a_i^2 n_i + A^2,
    x' A x  = sum_i a_i^2 (3 n_i - 1) + 2 A sum_i a_i
            = sum_i a_i^2 (3 n_i - 1) + 2 A^2,

and  <x,1> = 0  iff  sum_i a_i (n_i + 1) = 0.

So lambda_2(A(G)) >= R(n_1,n_2,n_3), the largest generalised eigenvalue of
(M, K) restricted to the 2‑dimensional space c^perp, where

    M = diag(3 n_i - 1) + 2 1 1',   K = diag(n_i) + 1 1',   c = (n_i + 1).

Sanity check: the two‑arm bottleneck vector (1/n_1) 1_{C_1} - (1/n_2) 1_{C_2}
lies in W, so R >= 3 - (n_1^2+n_2^2)/(n_1 n_2 (n_1+n_2)).
'''  # noqa: D400
from __future__ import annotations

import numpy as np
from typing import Tuple

RHO3: float = 2 * np.sqrt(2)


def _orthonormal_basis_orthogonal_to(v: np.ndarray) -> np.ndarray:
    """Return an orthonormal basis of the subspace orthogonal to *v*.

    Parameters
    ----------
    v:
        A 1‑D array of length 3.  The function normalises *v* and then uses a QR
        decomposition to obtain a 3×2 matrix whose columns form an orthonormal
        basis of ``v``\ :math:`^{\perp}`.

    Returns
    -------
    np.ndarray
        A 3×2 array ``B`` with ``B.T @ B == I_2`` and ``v @ B == 0``.
    """
    if v.shape != (3,):
        raise ValueError("Input vector must be one‑dimensional of length 3.")
    # Normalise to avoid overflow/underflow in the QR step.
    v_norm = v / np.linalg.norm(v)
    # Stack the normalised vector with the identity; QR yields an orthogonal
    # matrix whose first column spans ``v_norm``.
    q, _ = np.linalg.qr(np.column_stack([v_norm, np.eye(3)]), mode="reduced")
    # The remaining two columns form the desired basis.
    return q[:, 1:]


def best_quotient(n1: int, n2: int, n3: int) -> float:
    """Return the maximal Rayleigh quotient ``R(n1,n2,n3)``.

    The implementation follows the description in the module docstring but
    uses a numerically stable formulation based on solving a symmetric
    generalized eigenvalue problem ``M x = λ K x`` after projecting onto the
    subspace orthogonal to ``c = (n_i+1)``.

    Parameters
    ----------
    n1, n2, n3:
        Positive odd integers ``>= 5`` representing the sizes of the three
        arms.  A :class:`ValueError` is raised for invalid inputs.

    Returns
    -------
    float
        The largest eigenvalue of the projected problem, i.e. ``R(n1,n2,n3)``.
    """
    # Validate inputs early – this guards against division by zero later.
    for name, n in zip(("n1", "n2", "n3"), (n1, n2, n3)):
        if not isinstance(n, int) or n < 5 or n % 2 == 0:
            raise ValueError(f"{name} must be an odd integer >= 5, got {n!r}")

    nn = np.array([n1, n2, n3], dtype=float)
    # Matrices from the theoretical formulation.
    M = np.diag(3 * nn - 1.0) + 2.0 * np.ones((3, 3), dtype=float)
    K = np.diag(nn) + np.ones((3, 3), dtype=float)
    # Basis for the orthogonal complement of c = nn + 1.
    P = _orthonormal_basis_orthogonal_to(nn + 1.0)
    # Project the matrices onto the 2‑dimensional subspace.
    Mp = P.T @ M @ P
    Kp = P.T @ K @ P
    # Solve the symmetric generalized eigenvalue problem Kp^{-1} Mp.
    # Using ``np.linalg.solve`` avoids explicit matrix inversion and improves
    # numerical stability.
    S = np.linalg.solve(Kp, Mp)
    # Ensure symmetry before calling ``eigvalsh`` (round‑off can introduce a
    # tiny skew‑symmetric part).
    S_sym = (S + S.T) / 2.0
    eigenvalues = np.linalg.eigvalsh(S_sym)
    # The maximal eigenvalue corresponds to the Rayleigh quotient.
    return float(eigenvalues.max())


def bottleneck(n1: int, n2: int) -> float:
    """Two‑arm bottleneck lower bound.

    Parameters
    ----------
    n1, n2:
        Positive odd integers ``>= 5``.
    """
    if n1 <= 0 or n2 <= 0:
        raise ValueError("Arm sizes must be positive integers.")
    # Use ``float`` to guarantee Python 3 division semantics.
    return 3.0 - (float(n1 * n1 + n2 * n2) / (n1 * n2 * (n1 + n2)))


def main() -> None:
    """Simple CLI demonstration used in the original script.

    It prints a few sanity checks and then scans a range of odd triples up to
    300, reporting the minimal Rayleigh quotient encountered.
    """
    print(f"2*sqrt(2) = {RHO3:.9f}")
    b = bottleneck(5, 5)
    print(
        f"sanity: bottleneck(5,5) = {b:.9f}   R(5,5,5) = {best_quotient(5,5,5):.9f}"
        f"   (R must be >= bottleneck)"
    )
    print()
    print(f"{'n1':>5} {'n2':>5} {'n3':>6} {'bottleneck':>12} {'R':>14}   >2sqrt2?")
    cases = [
        (5, 5, 5), (5, 5, 7), (5, 5, 9), (5, 5, 11), (5, 5, 101), (5, 5, 9999),
        (5, 7, 7), (5, 7, 9), (5, 7, 9999), (5, 9, 9), (7, 7, 7), (7, 9, 9999),
        (5, 11, 9999), (5, 21, 9999), (5, 51, 9999), (5, 101, 9999),
        (5, 9999, 9999), (9999, 9999, 9999), (5, 5, 101), (5, 5, 10001),
    ]
    for n1, n2, n3 in cases:
        R = best_quotient(n1, n2, n3)
        print(
            f"{n1:5d} {n2:5d} {n3:6d} {bottleneck(n1, n2):12.9f} {R:14.9f}"
            f"   {'YES' if R > RHO3 else 'NO'}"
        )

    print("\nscan over odd 5 <= n1 <= n2 <= n3 <= 300:")
    odds = list(range(5, 301, 2))
    worst, at = float('inf'), None
    bad: list[Tuple[int, int, int, float]] = []
    for i, n1 in enumerate(odds):
        for j in range(i, len(odds)):
            n2 = odds[j]
            for k in range(j, len(odds)):
                n3 = odds[k]
                R = best_quotient(n1, n2, n3)
                if R < worst:
                    worst, at = R, (n1, n2, n3)
                if R <= RHO3 + 1e-12:
                    bad.append((n1, n2, n3, R))
    print(f"  min R = {worst:.9f} at {at}")
    print(f"  # scanned triples with R <= 2 sqrt 2: {len(bad)}")
    for x in bad[:20]:
        print("   ", x)


if __name__ == "__main__":
    main()
