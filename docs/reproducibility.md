# Reproducibility

## Requirements

| Component | Version used | Notes |
|---|---|---|
| Python | 3.14.6 | 3.10+ is sufficient |
| numpy | 2.4.6 | `numpy.linalg.eigvalsh` for all spectra |
| networkx | 3.6.1 | blossom algorithm for perfect matchings; `node_connectivity` |
| igraph | 1.0.0 | only in the (non-shipped-as-critical) cubic generator |
| sympy | 1.14.x | exact symbolic derivation of $M_\ell,K_\ell$ |
| OS | Windows (PowerShell) | nothing platform-specific; commands are plain `python` |

No compilation, no GPU, no network access, no proprietary software, no random
seeds (every computation is deterministic).

## Commands

```powershell
python tests/test_project.py              # 18 tests; expect "18/18 tests passed"
python experiments/run_all_experiments.py # Tables 1 and 2 of the paper
python experiments/case_b_symbolic.py      # symbolic derivation + cross-validation
python experiments/case_b_witness.py      # the (P,R) determinant and the l = 5 witness
python experiments/cut_vertex_two_case.py # Case A / Case B sanity checks
```

Total runtime: under one minute on the machine above.

## Expected output

`python tests/test_project.py` ends with

```
18/18 tests passed
```

`python experiments/run_all_experiments.py` prints

```
2*sqrt(2) = 2.828427124746
2*sqrt(d-1) for d=3; Marcus--Spielman--Srivastava bound on min_sigma lam_max(A_sigma) = 2.828427
Table 1: small cubic graphs
graph                        n        lam2  Ramanujan     beta(G)    lmax_min  kappa  girth
K_4                          4   -1.000000       True    2.236068    1.000000      3      3
K_{3,3}                      6    0.000000       True    2.000000    2.000000      3      3
triangular prism C_3 x K_2   6    1.000000       True    2.236068    2.000000      3      3
cube C_4 x K_2               8    1.000000       True    1.732051    1.732051      3      3
pentagonal prism            10    1.618034       True    2.236068    2.236068      3      3
hexagonal prism             12    2.000000       True    2.000000    2.000000      3      3
Heawood                     14    1.414214       True    2.449490    2.449490      3      5
Möbius--Kantor GP(8,3)      16    1.732051       True    2.414214    2.414214      3      5
```

`python experiments/case_b_symbolic.py` ends with

```
cross-validation: PASSED
```

## What is verified computationally and what is proved

**Proved by hand** (see `mathematics/proofs.md`):

* Lemma 0 (structure at a cut vertex);
* Theorem 1 (two-arm bottleneck bound) — pure Rayleigh-quotient algebra;
* Corollary 1 ($\lambda_2\ge20/7>2\sqrt2$ when two arms have order $\ge7$);
* Theorem 2 (cubic graphs with a cut vertex are not Ramanujan);
* Theorem 3, sub-case $\ell\ge7$: the sign of
  $-4(\ell+11)((14\sqrt2-19)\ell-5\sqrt2+2)$, which is exact algebra;
* Theorem 3, sub-case $\ell=5$: the exact inequality
  $19584^2-2\cdot13740^2=6\,017\,856>0$.

**Machine-checked finite verification** (these are *not* the proof, but they
guard it):

* the pencil $(M_\ell,K_\ell)$ of Theorem 3, re-derived symbolically from the
  three combinatorial identities and compared with matrices built from explicit
  coordinate vectors on concrete graphs, for $\ell=5,7,9,11,21,101,1001$
  (agreement to $10^{-9}$);
* $\lambda^*(\ell)>2\sqrt2$ for every odd $\ell\le199$ (double precision;
  the smallest margin is $2.285\times10^{-2}$ at $\ell=5$, far above
  double-precision error);
* the determinant sign for every odd $\ell\in[7,399]$;
* $\beta(G)$ for the eight graphs of Table 1, exact over switching classes.

## Determinism

No random number generation is used anywhere in the shipped code paths except
in `tests/test_project.py::test_switching_invariance`, which uses
`numpy.random.default_rng(0)` with a fixed seed.

## Known limitations of the code

* `src/signings/cubic.py` does **not** deduplicate correctly: it uses
  `igraph.Graph.canonical_permutation()`, which is not a canonical form (two
  isomorphic labelings can receive different permutations). It is retained only
  because it reproduces the classical counts of connected cubic graphs for
  $n=4,6,8,10$; it is used in no theorem and in no table. See
  `research/proof_attempts.md` §B3.
* `src/signings/core.py::vertex_connectivity` falls back to a brute-force
  routine if networkx is unavailable; the fallback is slow but was not used in
  any reported number.
* All reported spectra are in double precision. Margins to the decision
  thresholds ($\sqrt2$, $2\sqrt2$, $20/7$) are at least $10^{-2}$ everywhere, so
  rounding cannot change any conclusion.