# Known results: what is proved here and what is borrowed

## PROVED HERE

| # | Statement | File |
|---|---|---|
| Lemma 0 | Structure of a connected cubic graph at a cut vertex: $G-z$ has exactly 3 components, each odd of order $\ge5$, each attached to $z$ by one edge | `mathematics/proofs.md` |
| Theorem 1 | Two-arm bottleneck bound $\lambda_2(A(G))\ge\max_{i<j}\left(3-\frac{n_i^2+n_j^2}{n_in_j(n_i+n_j)}\right)$ | `mathematics/proofs.md` |
| Corollary 1 | If two of the three arm orders are $\ge7$, then $\lambda_2(A(G))\ge\frac{20}{7}>2\sqrt2$ | `mathematics/proofs.md` |
| Theorem 2 | Every connected cubic graph with a cut vertex satisfies $\lambda_2(A(G))>2\sqrt2$, hence is not Ramanujan | `mathematics/proofs.md` |
| Corollary 2 | Every connected Ramanujan cubic graph is 2-connected | `mathematics/proofs.md` |
| Theorem 3 | For arm orders $(5,5,\ell)$, the 4-dimensional space $\operatorname{span}\{p_1+p_2,q_1+q_2,p_3,q_3\}$ yields $\lambda_2\ge\lambda^*(\ell)>2\sqrt2$ for all odd $\ell\ge5$; explicit pencil given | `mathematics/proofs.md` |
| Proposition 4 | If a connected Ramanujan $d$-regular graph $G$ is a 2-lift of a $d$-regular graph $G_0$ via $\sigma_0$, then $G_0$ is connected, $\mathrm{Spec}(A_{\sigma_0}(G_0))\subseteq[-\rho_d,\rho_d]\cup\{-d\}$, and $-d$ cannot occur when $G_0$ is bipartite | `mathematics/proofs.md` |

## ESTABLISHED IN PRIOR WORK (verified citations)

| Statement | Source |
|---|---|
| Bilu--Linial signing conjecture $\forall$ $d$-regular $G\;\exists\sigma:\|A_\sigma\|\le2\sqrt{d-1}$ — **false** | Z. Xu, arXiv:2609.15591v4 (2026) |
| Every $d$-regular graph admits $\sigma$ with all signed adjacency eigenvalues $\le2\sqrt{d-1}$ (one-sided) | Marcus--Spielman--Srivastava, *Interlacing families I*, Ann. of Math. 182 (2015) Thm 5.3 |
| Every $d$-regular graph admits $\sigma$ with $\|A_\sigma\|<2\sqrt{2(d-1)}$ | Ravichandran--Srivastava, IMRN 2021 |
| Every $d$-regular graph ($d\ge2$) admits $\sigma$ with $\|A_\sigma\|<\frac{3+\sqrt5}{2}\sqrt{d-1}$ | Lin--Zhou, arXiv:2609.15715 (2026) |
| Every connected $d$-regular graph ($d\ge3$) has an $\ell$-lift with new eigenvalues $\le2\sqrt{d-1}$ | Hall--Puder--Sawin, Adv. Math. 323 (2018) |
| Bilu--Linial conjecture for bipartite $d$-regular graphs | follows from MSS (spectrum of $A_\sigma$ symmetric for bipartite $G$) |

## COMPUTATIONALLY VERIFIED (finite domains, exact or double precision, as stated)

| Claim | Domain | Script |
|---|---|---|
| $\beta(K_4)=\sqrt5$, $\beta(K_{3,3})=2$, $\beta(\text{cube})=\sqrt3$ | exhaustive over $2^{m-n+1}$ switching classes | `experiments/run_all_experiments.py` |
| $\beta(G)\le3$ for the eight cubic graphs of Table 1 | same | same |
| The pencil of Theorem 3 computed symbolically from the three combinatorial identities agrees with matrices built from explicit coordinate vectors | $\ell=5,7,9,11,21,101,1001$; agreement to $10^{-9}$ | `tests/test_project.py` |
| $\lambda^*(l)>2\sqrt2$ for all odd $\ell\le199$ | odd $\ell\in[5,199]$ | `tests/test_project.py` |
| The $(P,R)$ principal minor of $M_\ell-2\sqrt2K_\ell$ is negative for all odd $\ell\ge7$, $\ell<400$ | odd $\ell\in[7,399]$ | `tests/test_project.py` |
| No three-arm graph with arm orders in the tested list is Ramanujan | 10 arm-order triples, exact eigenvalues | `tests/test_project.py` |
| Counts of connected cubic graphs generated match A002851 | $n=4,6,8,10$ | `experiments/_gen_check.py` (deleted; see `research/proof_attempts.md` B3) |
| Theorem 2 holds on random three-arm cubic graphs | 4000 random graphs, arm orders $5..11$, seed 20261013; 0 Ramanujan, smallest margin $\lambda_2-2\sqrt2=0.0273$ | `experiments/verify_theorem2_random.py` |
| Proposition 4 holds | all $37952$ signings of $K_4$, $K_{3,3}$, the triangular prism, the cube, the Petersen graph | `experiments/verify_lift_proposition.py` |

## WITHDRAWN CLAIMS (recorded so they are not mistaken for results)

| Claim | Why withdrawn |
|---|---|
| "Vertex sums are invariant under vertex switching" | false; see `research/proof_attempts.md` B1 |
| "$\max_v\|\sigma(v)\|$ bounds $\|A_\sigma\|$" | false: $K_4$ with a perfect-matching signing | `research/proof_attempts.md` A1/B4 |
| "$\beta(G_0)\le\rho_d$ for every base of a Ramanujan 2-lift" | false: 552 counterexample lifts | `research/proof_attempts.md` B5 |
| "the exceptional signed eigenvalue occurs iff the base is non-bipartite" | false: 32240 counterexamples | `research/proof_attempts.md` B5 |

## CONJECTURED (not proved here)

* Conjecture R — Ramanujan signing for regular Ramanujan graphs. **Open.**
* Conjecture R3 — the $d=3$ case. **Open.**
* (Side observation, CONJECTURED.) Every *3-connected* cubic graph satisfies
  $\beta(G)\le2\sqrt2$.  This is suggested by (a) Theorem 2 and (b) the absence
  of any cubic graph with $\beta>2\sqrt2$ among the small graphs examined, but
  no proof is given and no exhaustive census beyond $n\le16$ was run here.

## UNRESOLVED

* Whether Ramanujan cubic graphs admit Ramanujan signings.
* The smallest $d$ for which $2\sqrt{d-1}$ fails, and the smallest cubic graph
  with $\beta>2\sqrt2$.
* Whether the $(5,5,\ell)$ obstruction found in Theorem 3 is the *only* type of
  obstruction that could defeat Ramanujan signing.