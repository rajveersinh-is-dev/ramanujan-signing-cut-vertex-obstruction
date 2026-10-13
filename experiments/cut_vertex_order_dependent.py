"""Order-dependent uniform test subspace for cubic graphs with a cut vertex.

W(m,k,l) = span{ 1_{C_i} - n_i e_z ,  e_{t_i} - (1/n_i) 1_{C_i} : i = 1,2,3 }.

All six generators are orthogonal to the all-ones vector, and the subspace
depends on the arm orders (n_1,n_2,n_3) only -- not on the internal structure of
the arms.  Hence

    lambda_2(A(G))  >=  max_{x in W, x != 0}  x'A(G)x / ||x||^2

is a genuine uniform bound over all cubic graphs with a cut vertex whose arm
orders are (m, k, l).
"""
from __future__ import annotations


import numpy as np


RHO3 = 2 * np.sqrt(2)


def arm_of_order(n) -> tuple:
    """Arm of order.
    
    Args:
        n:
    
    Returns:
        tuple: Result of type tuple
    
    """
    assert n % 2 == 1 and n >= 5
    if n == 5:
        he = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    else:
        kk = (n - 1) // 2
        he = sorted(set([(i, (i + 1) % kk) for i in range(kk)]
                        + [(i + kk, ((i + 1) % kk) + kk) for i in range(kk)]
                        + [(i, i + kk) for i in range(kk)]))
    E = [e for e in he if frozenset(e) != frozenset((0, 1))]
    E = E + [(n - 1, 0), (n - 1, 1)]
    return n, sorted(E), n - 1


def build(orders) -> tuple:
    """Create the object.
    
    Args:
        orders:
    
    Returns:
        tuple: Result of type tuple
    
    """
    off = 0
    edges, tips = [], []
    for o in orders:
        n, e, t = arm_of_order(o)
        edges += [(u + off, v + off) for u, v in e]
        tips.append(off + t)
        off += n
    z = off
    edges += [(t, z) for t in tips]
    n = z + 1
    A = np.zeros((n, n))
    for u, v in edges:
        A[u, v] = A[v, u] = 1
    offs, o2 = [], 0
    for o in orders:
        offs.append(o2)
        o2 += o
    return n, A, offs, z, tips


def uniform_bound(orders, exact_check=False) -> tuple:
    """Uniform bound.
    
    Args:
        orders:
        exact_check (bool):
    
    Returns:
        tuple: Result of type tuple
    
    """
    n, A, offs, z, tips = build(orders)
    cols = []
    for i, o in enumerate(orders):
        v = np.zeros(n)
        v[offs[i]:offs[i] + o] = 1
        v[z] = -o
        cols.append(v)
        w = np.zeros(n)
        w[tips[i]] += 1
        w[offs[i]:offs[i] + o] -= 1.0 / o
        cols.append(w)
    V = np.column_stack(cols)
    U, s, _ = np.linalg.svd(V, full_matrices=False)
    Q = U[:, :int((s > 1e-9).sum())]
    assert abs(Q.sum(axis=0)).max() < 1e-8
    B = Q.T @ (A @ Q)
    val = float(np.linalg.eigvalsh((B + B.T) / 2).max())
    if exact_check:
        return val, float(np.linalg.eigvalsh(A)[::-1][1])
    return val


def main():
    """Entry point — parse arguments and run the main computation.
    
    """
    print(f"2*sqrt(2) = {RHO3:.12f}\n")
    print("the hard family: arms (5,5,l)")
    for l in (5, 7, 9, 11, 13, 15, 17, 19, 21, 25, 31, 51, 101, 201, 501,
              1001, 5001, 10001, 100001):
        b, lam2 = uniform_bound((5, 5, l), exact_check=True)
        print(f"  l={l:7d}  bound={b:.9f}  lam2={lam2:.9f}  "
              f"{'YES' if b > RHO3 else 'NO'}")
    print()
    print("scan over odd 5 <= n1 <= n2 <= n3 <= 120:")
    odds = list(range(5, 121, 2))
    worst, at = 1e18, None
    bad = []
    for i, n1 in enumerate(odds):
        for j in range(i, len(odds)):
            for k in range(j, len(odds)):
                b = uniform_bound((n1, odds[j], odds[k]))
                if b < worst:
                    worst, at = b, (n1, odds[j], odds[k])
                if b <= RHO3 + 1e-12:
                    bad.append((n1, odds[j], odds[k], b))
    print(f"  min bound = {worst:.12f} at {at}")
    print(f"  #triples with bound <= 2sqrt2: {len(bad)}")
    for x in bad[:
        20]:
        print("   ", x)


if __name__ == "__main__":
    main()