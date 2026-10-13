"""Finding a uniform exclusion subspace for the residual case n_1 = n_2 = 5.

We look for a small, order-dependent subspace of {x : <x,1> = 0} contained in
the graph on {z} U C_1 U C_2 (plus optionally C_3) whose top Rayleigh quotient
exceeds 2 sqrt 2 for every third-arm order l >= 5.

Generators: p_i = 1_{C_i} - n_i e_z,  q_i = e_{t_i} - (1/n_i) 1_{C_i},
and the "symmetrised" versions p+ = p_1+p_2, q+ = q_1+q_2.
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
    return n, sorted(E + [(n - 1, 0), (n - 1, 1)]), n - 1


def build(orders) -> tuple:
    """Create the object.
    
    Args:
        orders:
    
    Returns:
        tuple: Result of type tuple
    
    """
    off, edges, tips = 0, [], []
    for o in orders:
        nn, e, t = arm_of_order(o)
        edges += [(u + off, v + off) for u, v in e]
        tips.append(off + t)
        off += nn
    z = off
    edges += [(t, z) for t in tips]
    N = z + 1
    A = np.zeros((N, N))
    for u, v in edges:
        A[u, v] = A[v, u] = 1
    offs, o2 = [], 0
    for o in orders:
        offs.append(o2)
        o2 += o
    return N, A, offs, z, tips


def gens(N, offs, z, tips, orders) -> tuple:
    """Gens.
    
    Args:
        N:
        offs:
        z:
        tips:
        orders:
    
    Returns:
        tuple: Result of type tuple
    
    """
    P, Q = [], []
    for i, o in enumerate(orders):
        p = np.zeros(N); p[offs[i]:offs[i] + o] = 1; p[z] = -o
        q = np.zeros(N); q[tips[i]] += 1; q[offs[i]:offs[i] + o] -= 1.0 / o
        P.append(p); Q.append(q)
    return P, Q


def top_of_span(A, vectors) -> float:
    """Top of span.
    
    Args:
        A:
        vectors:
    
    Returns:
        float: Result of type float
    
    """
    V = np.column_stack(vectors)
    assert abs(V.sum(axis=0)).max() < 1e-8
    K = V.T @ V
    M = V.T @ (A @ V)
    ev, U = np.linalg.eigh(K)
    Ki = U @ np.diag(1.0 / np.sqrt(ev)) @ U.T
    S = (Ki @ M @ Ki)
    return float(np.linalg.eigvalsh((S + S.T) / 2).max())


def pencil(A, vectors) -> tuple:
    """Pencil.
    
    Args:
        A:
        vectors:
    
    Returns:
        tuple: Result of type tuple
    
    """
    V = np.column_stack(vectors)
    K = V.T @ V
    M = V.T @ (A @ V)
    return M, K


def main():
    """Entry point — parse arguments and run the main computation.
    
    """
    ls = (5, 7, 9, 11, 15, 21, 25, 31, 51, 101, 501, 5001)
    schemes = {
        "p+ , q+":                       lambda P, Q: [P[0] + P[1], Q[0] + Q[1]],
        "p+ , q+ , p3":                  lambda P, Q: [P[0] + P[1], Q[0] + Q[1], P[2]],
        "p+ , q+ , q3":                  lambda P, Q: [P[0] + P[1], Q[0] + Q[1], Q[2]],
        "p+ , q+ , p3 , q3":             lambda P, Q: [P[0] + P[1], Q[0] + Q[1], P[2], Q[2]],
        "p+ , q1-q2":                    lambda P, Q: [P[0] + P[1], Q[0] - Q[1]],
        "p+ , q+ , p1-p2":               lambda P, Q: [P[0] + P[1], Q[0] + Q[1], P[0] - P[1]],
        "p1,p2,p3":                      lambda P, Q: [P[0], P[1], P[2]],
        "p+,q+,p3,q3,q1-q2":             lambda P, Q: [P[0] + P[1], Q[0] + Q[1], P[2], Q[2], Q[0] - Q[1]],
        "full 6":                        lambda P, Q: [P[i] for i in range(3)] + [Q[i] for i in range(3)],
    }
    print(f"2*sqrt(2) = {RHO3:.9f}\n")
    print(f"{'scheme':<24}" + "".join(f"{('l=' + str(l)):>12}" for l in ls[:8]) + "   min")
    for name, sch in schemes.items():
        vals = []
        for l in ls:
            N, A, offs, z, tips = build((5, 5, l))
            P, Q = gens(N, offs, z, tips, (5, 5, l))
            vals.append(top_of_span(A, sch(P, Q)))
        print(f"{name:<24}" + "".join(f"{v:12.6f}" for v in vals[:8])
              + f"   {min(vals):.6f}")

    print("\n2x2 pencil for span{p+ , q+} on the 16-vertex (5,5,5) graph:")
    N, A, offs, z, tips = build((5, 5, 5))
    P, Q = gens(N, offs, z, tips, (5, 5, 5))
    M, K = pencil(A, [P[0] + P[1], Q[0] + Q[1]])
    print("M =\n", np.round(M, 6))
    print("K =\n", np.round(K, 6))

    print("\n4x4 pencil for span{p+, q+, p3, q3} on (5,5,5):")
    M, K = pencil(A, [P[0] + P[1], Q[0] + Q[1], P[2], Q[2]])
    print("M =\n", np.round(M, 6))
    print("K =\n", np.round(K, 6))


if __name__ == "__main__":
    main()