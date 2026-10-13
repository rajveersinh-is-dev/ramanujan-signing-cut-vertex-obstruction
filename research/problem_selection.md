# Problem selection

## The problem attacked

**Ramanujan signing of regular graphs.**  For a finite connected simple $d$-regular
graph $G$ and a signing $\sigma\colon E(G)\to\{\pm1\}$ write $A_\sigma(G)$ for the
signed adjacency matrix and $\beta(G)=\min_\sigma\|A_\sigma(G)\|$.  Let
$\rho_d=2\sqrt{d-1}$ (the Ramanujan bound for degree $d$).

> **Conjecture R (stated as open by Z. Xu, arXiv:2609.15591v4, 29 Sep 2026).**
> If $G$ is a connected $d$-regular *Ramanujan* graph then $\beta(G)\le\rho_d$.

The present project attacks its $d=3$ instance

> **Conjecture R3.**  If $G$ is a connected cubic graph with
> $\lambda_2(A(G))\le2\sqrt2$ then $\beta(G)\le2\sqrt2$,

and establishes the structural results 1–3 in `mathematics/main_results.md`.

## Why this problem

**It is open and its frontier is very recent.**  The original Bilu–Linial
signing conjecture (2006) was *disproved* on 14 Sep 2026 by Z. Xu
(arXiv:2609.15591v4), and the same paper's conclusion lists exactly two open
questions, the first of which is Conjecture R.  At the time of writing (Oct 2026)
this is a matter of weeks, not decades — a good sign that the gap is still
genuinely open rather than already settled by later work.

**The decisive work is mathematical, not computational.**  Signing norms of
small graphs are computable in a few thousand operations, but the *statement*
concerns all graphs of arbitrary order.  Nothing in the mathematics below needs
a large machine: the results proved here are elementary, self-contained and
consist of test-vector constructions plus exact arithmetic.

**There is a clean, checkable "structural" target.**  Xu's counterexample is
built around a cut vertex (his core $J$ has a central vertex whose removal
leaves three huge components).  It is therefore natural to ask whether
cut-vertex structure can ever obstruct *Ramanujan* signings.  Theorem 2 answers
this negatively: cut vertices in cubic graphs automatically destroy Ramanujan-ness.
This is a genuine reduction of where a counterexample to Conjecture R3 could
possibly live.

**Verification is cheap and exact.**  A cubic graph on $n$ vertices has
$2^{n/2+1}$ switching classes, so $\beta(G)$ is exactly computable for $n\lesssim20$
by exhaustive search; all the numerical claims in this repository are produced
in seconds.

## Alternatives considered and rejected

| Candidate | Why rejected |
|---|---|
| The Bilu–Linial conjecture $\beta(G)\le2\sqrt{d-1}$ for all $d$-regular $G$ | **Disproved** (Xu, 29 Sep 2026). |
| Erdos--Nathanson: minimal asymptotic bases of order $h$ have density 0 | Open but the known partial results require either strong structural hypotheses or analytic number theory; no plausible short attack. |
| Krausz--Schmagerl "little cross number conjecture" | Open, but proofs in the literature are 30+ pages of intricate unit-equation combinatorics. |
| Zero-sum Ramsey problem $Z(v)=p$ (Gao--Geroldinger) | The decisive lemmas are arithmetic (unit equations in small numbers of variables); not a good fit for structural/combinatorial methods. |
| Bilu--Linial restricted to Ramanujan base graphs (Q1 of Xu) | This *is* Conjecture R; attacked here. |
| Lin--Zhou improved universal bound $\beta\le\frac{3+\sqrt5}{2}\sqrt{d-1}$ | Would require redoing interlacing-family machinery; out of scope. |
| Finding a small cubic graph with $\beta>2\sqrt2$ | Genuine open question, but the answer is a minimum-over-a-census computation; would make the project computational rather than mathematical. Recorded as an open direction instead. |

## Outcome

The full conjectures (R and R3) are **not** solved.  What is established here is
Theorem 1 (two-arm bottleneck bound), Corollary 1 ($\lambda_2\ge20/7>2\sqrt2$ when
two arms have order $\ge7$), Theorem 3 (a uniform four-dimensional test space for
the residual case) and Theorem 2 (a connected cubic graph with a cut vertex is
never Ramanujan).  See `README.md` for the exact statement and
`research/proof_attempts.md` for the approaches that failed.