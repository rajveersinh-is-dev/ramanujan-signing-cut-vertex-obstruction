# Ramanujan signing of regular graphs: a structural obstruction and the state of the open problem

A research project on the **Ramanujan signing conjecture** (Xu, 2026), with a
new structural theorem for cubic graphs.

---

## The problem

For a finite connected simple $d$-regular graph $G$, a **signing** is a map
$\sigma\colon E(G)\to\{\pm1\}$ and $A_\sigma(G)$ is the signed adjacency matrix.
Write
$$
\beta(G):=\min_{\sigma}\ \|A_\sigma(G)\|_{\mathrm{op}},\qquad \rho_d:=2\sqrt{d-1}.
$$

Bilu and Linial (2006) conjectured $\beta(G)\le\rho_d$ for every $d$-regular
graph. **That conjecture is false**: Z. Xu constructed a cubic graph with
$\beta(G)>2\sqrt2$ (arXiv:2609.15591v4, 29 Sep 2026). Xu's paper ends by
stating, as an open question, the following.

> **Conjecture R (open).** *Ramanujan signing.* If $G$ is a connected $d$-regular
> **Ramanujan** graph — every non-trivial eigenvalue of $A(G)$ has modulus at
> most $\rho_d$ — then $\beta(G)\le\rho_d$.

This project attacks the $d=3$ case,

> **Conjecture R3 (open).** If $G$ is connected and cubic with
> $\lambda_2(A(G))\le2\sqrt2$, then $\beta(G)\le2\sqrt2$,

and proves the structural results below. **Neither conjecture is solved here.**

---

## What is proved here

All statements are proved in [`mathematics/proofs.md`](mathematics/proofs.md)
(mirrored in [`paper/main.tex`](paper/main.tex)).

> **Theorem 1 (two-arm bottleneck bound).** *PROVED HERE.* Let $G$ be connected
> and cubic, $z$ a cut vertex, and $G-z=C_1\dot\cup C_2\dot\cup C_3$ with
> $n_i=|C_i|$, $n_1\le n_2\le n_3$. Then
> $$\lambda_2(A(G))\ \ge\ \max_{i<j}\Bigl(3-\frac{n_i^2+n_j^2}{n_in_j(n_i+n_j)}\Bigr).$$

> **Corollary 1.** *PROVED HERE.* If at least two of $n_1,n_2,n_3$ are $\ge7$,
> then $\lambda_2(A(G))\ge\frac{20}{7}>2\sqrt2$.

> **Theorem 3 (residual case $(5,5,\ell)$).** *PROVED HERE.* If $n_1=n_2=5$ and
> $n_3=\ell$ (odd, $\ge5$), the four-dimensional space
> $\operatorname{span}\{\mathbf 1_{C_1}+\mathbf 1_{C_2}-10e_z,\;
> (e_{t_1}-\tfrac15\mathbf 1_{C_1})+(e_{t_2}-\tfrac15\mathbf 1_{C_2}),\;
> \mathbf 1_{C_3}-\ell e_z,\; e_{t_3}-\tfrac1\ell\mathbf 1_{C_3}\}$
> gives $\lambda_2(A(G))\ge\lambda^*(\ell)>2\sqrt2$ for every odd $\ell\ge5$,
> with an explicit pencil $(M_\ell,K_\ell)$ depending only on the arm orders.
> Numerically $\lambda^*(5)=2.8508$ and $\lambda^*(\ell)\nearrow2.9561$.

> **Theorem 2 (main result).** *PROVED HERE.* **Every connected cubic graph with
> a cut vertex is not Ramanujan:**
> $$\lambda_2(A(G))>2\sqrt2 .$$

> **Corollary 2.** *PROVED HERE.* Every connected Ramanujan cubic graph is
> 2-connected.

> **Proposition 4 (2-lifts).** *PROVED HERE, after one correction.* If a
> connected Ramanujan $d$-regular graph $G$ is a 2-lift of a $d$-regular graph
> $G_0$ via the signing $\sigma_0$, then $G_0$ is connected,
> $\mathrm{Spec}(A_{\sigma_0}(G_0))\subseteq[-\rho_d,\rho_d]\cup\{-d\}$, and the
> exceptional value $-d$ cannot occur when $G_0$ is bipartite. In particular a
> *lift-compatible* signing of a 2-lift is never Ramanujan (the base spectrum
> always contains $d$), while a bipartite base always does have a Ramanujan
> signing.

### Why this is worth knowing

Xu's counterexample to Bilu--Linial is built around a cut vertex: his core graph
$J$ has a central vertex whose removal leaves three enormous components.
Theorem 2 shows that this mechanism can **never** be used to refute Conjecture
R for $d=3$: any cubic graph with such a bottleneck automatically fails the
Ramanujan hypothesis. So the obstruction is genuinely two-sided — the graphs
that break the signing bound are not Ramanujan, and the graphs that are Ramanujan
cannot have this kind of bottleneck.

The proof is entirely structural: a Rayleigh-quotient test vector (the classic
"bottleneck" construction) in the easy case, and a four-dimensional test space
reduced by an exact $4\times4$ pencil analysis in the one remaining
configuration. **No large computation is involved anywhere in the proof.**

---

## Open status

| Statement | Status | Evidence |
|---|---|---|
| Bilu--Linial signing conjecture | **FALSE** (disproved 2026-09-14, v4 2026-09-29) | Z. Xu, *A 3-regular counterexample to the Bilu--Linial signing conjecture*, arXiv:2609.15591v4; abstract and §8 read on 2026-10-13 |
| MSS one-sided bound $\min_\sigma\lambda_{\max}(A_\sigma)\le2\sqrt{d-1}$ | **ESTABLISHED IN PRIOR WORK** | Marcus–Spielman–Srivastava, *Interlacing families I*, Ann. of Math. 182 (2015) 307–325, Thm 5.3 |
| Universal norm bound $\beta(G)<\tfrac{3+\sqrt5}{2}\sqrt{d-1}$ | **ESTABLISHED IN PRIOR WORK** | Lin–Zhou, arXiv:2609.15715 (2026) |
| **Conjecture R (Ramanujan signing)** | **OPEN** | stated as open in Xu (2026) §8; no later resolution found in arXiv searches up to 2026-10-13 |
| **Conjecture R3 ($d=3$)** | **OPEN** | as above |
| Every *3-connected* cubic graph has $\beta\le2\sqrt2$ | **CONJECTURED** | suggested by Theorem 2; no proof here |

**Novelty assessment (conservative).** The two-arm Rayleigh estimate of
Theorem 1 is an elementary instance of the standard "bottleneck" / Cheeger-type
bound and should be regarded as **classical technique applied carefully**, not
as a new discovery. Theorem 3's pencil and the case analysis of Theorem 2 were
derived here and I did not find an equivalent statement in the literature
searched, but I classify the result as **potentially new, novelty uncertain**
rather than new: I cannot exclude that it appears as folklore. See
`research/literature_review.md` for exactly what was searched.

---

## Repository map

```
README.md                     this file
LICENSE                       MIT
CHANGELOG.md                  chronological record of commits
mathematics/
  definitions.md              signings, switching, beta(G), Ramanujan condition
  main_results.md             Theorem 1, Corollary 1, Theorem 2, Theorem 3, Conjecture R
  proofs.md                   complete audited proofs
paper/
  main.tex                    the paper (LaTeX)
  references.bib              bibliography
research/
  problem_selection.md        what was chosen, what was rejected, why
  literature_review.md        verified references and status evidence
  known_results.md            evidence labels for every claim
  proof_attempts.md           dead ends and CORRECTED errors
  research_log.md             chronological log, including corrections
src/signings/
  core.py                     signings, switching classes, spectra, exact beta(G)
  cubic.py                    exhaustive generator for connected cubic graphs
  theorems.py                 executable transcriptions of the proof steps
tests/test_project.py         18 tests
experiments/                  exploratory scripts; run_all_experiments.py is the
                              one that reproduces the numbers in the paper
docs/reproducibility.md       exact commands, versions, expected output
```

---

## Reproducing the results

Requirements: Python $\ge3.10$ with `numpy`, `networkx`, `igraph`, `sympy`.
No compilation, no GPU, no network access needed. Tested with Python 3.14.6,
numpy 2.4.6, networkx 3.6.1, igraph 1.0.0, sympy (1.14 era) on Windows.

```powershell
# the test suite (also re-derives the pencil of Theorem 3 and cross-checks it)
python tests/test_project.py

# the tables reported in the paper
python experiments/run_all_experiments.py

# the symbolic derivation of M_l, K_l from the three combinatorial identities
python experiments/case_b_symbolic.py

# independent randomised check of Theorem 2 (4000 random three-arm graphs)
python experiments/verify_theorem2_random.py

# exhaustive check of Proposition 4 (all 37952 signings of five base graphs)
python experiments/verify_lift_proposition.py
```

Expected: `18/18 tests passed`, the two tables of
`experiments/run_all_experiments.py`, `cross-validation: PASSED`,
`RESULT: PASS` from both verification scripts. Total runtime: under two
minutes.

---

## Computations actually run (and their scope)

Every number quoted as "computed" in this repository comes from
`experiments/run_all_experiments.py` (IEEE-754 double precision via
`numpy.linalg.eigvalsh`) or from the exact-rational/sympy scripts.

```
2*sqrt(2) = 2.828427124746

Table 1: small cubic graphs
graph                        n        lam2  Ramanujan     beta(G)    lmax_min  kappa  girth
K_4                          4   -1.000000       True    2.236068    1.000000      3      3
K_{3,3}                      6    0.000000       True    2.000000    2.000000      3      3
triangular prism C_3 x K_2   6    1.000000       True    2.236068    2.000000      3      3
cube C_4 x K_2               8    1.000000       True    1.732051    1.732051      3      3
pentagonal prism            10    1.618034       True    2.236068    2.236068      3      3
hexagonal prism             12    2.000000       True    2.000000    2.000000      3      3
Heawood                     14    1.414214       True    2.449490    2.449490      3      5
Mobius--Kantor GP(8,3)       16    1.732051       True    2.414214    2.414214      3      5
```

(`kappa` is the vertex connectivity; `lmax_min` is
$\min_\sigma\lambda_{\max}(A_\sigma)$, which MSS bounds by $2\sqrt{d-1}=2.828427$
— satisfied with room to spare in every row.)

Scope and limitations:

* $\beta(G)$ is **exact**, not sampled: for a connected cubic graph on $n$
  vertices there are exactly $2^{n/2+1}$ vertex-switching classes, all of which
  are enumerated. The listed graphs have $n\le16$, i.e. at most $512$ classes.
* All eight listed graphs are Ramanujan and all satisfy $\beta\le2\sqrt2$, as
  Conjecture R3 predicts. This is evidence on eight graphs, not a proof.
* $\lambda_2$ and the $A_\sigma$ eigenvalues are computed in double precision;
  the margins to $2\sqrt2$ here ($\ge0.4$) are far beyond floating-point error.
* The only *exact* arithmetic used in the proof itself is the evaluation of
  $19584^2-2\cdot13740^2=6\,017\,856>0$ and the symbolic sign of
  $-4(\ell+11)((14\sqrt2-19)\ell-5\sqrt2+2)$.

---

## Known limitations

* **The conjectures are open.** This is a partial result (Category C).
* Theorem 2 rules out cut vertices but says nothing about 2-vertex cuts,
  bridges that are not cut vertices of order 2, or bottlenecks that are
  "large but finite".
* The $d\ge4$ case of Conjecture R is untouched.
* The cubic graph generator in `src/signings/cubic.py` produces duplicates
  (its isomorphism filter uses `igraph.canonical_permutation()`, which is not a
  canonical form — see `research/proof_attempts.md` B3). It is validated for
  $n=4,6,8,10$ only up to the *number of distinct generated graphs*, and is not
  used in any theorem or table here.

## Remaining open questions

1. Conjecture R3: do Ramanujan cubic graphs admit signings with
   $\|A_\sigma\|\le2\sqrt2$?
2. Is every *3-connected* cubic graph 2-realisable, i.e. does it satisfy
   $\beta(G)\le2\sqrt2$?
3. What is the smallest cubic graph with $\beta(G)>2\sqrt2$? (Xu's example has
   astronomically many vertices; exhaustive search to $n=16$ is trivial with the
   switching-class machinery already in this repository.)
4. Does the bound of Theorem 1 have an equality case, and if so when?

---

## License

MIT — see `LICENSE`. Prior work cited in `paper/references.bib` remains under
its own licences; no prior text is reproduced here.