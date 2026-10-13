# Superseded exploratory scripts

These scripts are kept **as a record of mistakes made during the
investigation**.  They are not used by any theorem, table, or test, and at
least one of them is deliberately wrong.  Do not import from them.

| Script | Why superseded | See |
|---|---|---|
| `explore_cubic_signings.py` | First version of the signing code.  Its `signed_matrix` multiplied a vertex-indexed vector by an edge-indexed one and raised a `ValueError`; it also asserted that vertex sums are switching-invariant, which is false.  Replaced by `src/signings/core.py`. | `research/proof_attempts.md` B1 |
| `cut_vertex_pencil.py` | Contains the **hand-derived closed-form pencil** $(M_\ell,K_\ell)$ that is *wrong*: it disagrees with the pencil built from explicit coordinate vectors (the script prints `MISMATCH`).  It is retained because the mismatch is exactly how error B2 was detected — the wrong formula produces a top eigenvalue larger than $3$, impossible for a cubic graph.  The correct derivation is in `experiments/case_b_symbolic.py`. | `research/proof_attempts.md` B2 |

## Scripts that are current and run cleanly

```
arm_indicator_span.py        exploration: optimum on span of arm indicators
arm_mu.py                    exploration: mean-zero (second-order) test vectors
case_b_symbolic.py           CORRECT symbolic derivation of M_l, K_l + cross-validation
case_b_witness.py            the (P,R) determinant sign and the l = 5 rational witness
cut_vertex_enriched.py       exploration: test space enriched with mean-zero directions
cut_vertex_order_dependent.py  exploration: the 6-dimensional order-dependent space
cut_vertex_spectrum.py       exact lambda_2 for three-arm graphs, many arm triples
cut_vertex_two_case.py       Case A / Case B sanity checks (Case A proves (5))
cut_vertex_uniform.py        scan of the four-dimensional space over arm orders
find_case_b_subspace.py      which small test spaces work for the (5,5,l) family
residual_family.py           bound and top eigenvector across the (5,5,l) family
run_all_experiments.py       THE script that reproduces the paper's tables
```

`case_b_symbolic.py`, `case_b_witness.py`, `cut_vertex_two_case.py` and
`run_all_experiments.py` are the ones cited in `docs/reproducibility.md`.