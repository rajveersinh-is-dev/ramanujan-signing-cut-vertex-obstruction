# Changelog

All notable changes to this project, in order.

## [0.1.0] - 2026-10-13

### Research

- **research/problem_selection.md** — survey of candidate open problems;
  selection of the Ramanujan signing conjecture (Xu, 2026) and its cubic case.
- **research/literature_review.md** — verified status evidence for the Bilu–Linial
  signing conjecture (false as of 2026-09-29) and for Ramanujan signing (open);
  record of all searches performed.
- **research/known_results.md** — evidence labels for every claim in the project.
- **research/proof_attempts.md** — dead ends, and four errors made and corrected
  (switching-invariance of vertex sums; the $\max_v|\sigma(v)|$ bound; three
  faulty hand-derived bilinear forms; misuse of `igraph.canonical_permutation`).
- **research/research_log.md** — chronological log.

### Mathematics

- **mathematics/definitions.md** — signings, vertex switching, the switching-class
  count $2^{m-n+1}$, $\beta(G)$, the Ramanujan condition, the cut-vertex setup.
- **mathematics/main_results.md** — Theorem 1, Corollary 1, Theorem 2,
  Corollary 2, Theorem 3, and the two open conjectures.
- **mathematics/proofs.md** — complete, audited proofs of all of the above.

### Paper

- **paper/main.tex**, **paper/references.bib** — LaTeX write-up of the same
  results, with abstract, discussion, limitations and open questions.

### Code

- **src/signings/core.py** — signed adjacency matrices, switching classes,
  exact $\beta$, spectra, girth, vertex connectivity, perfect matchings.
- **src/signings/cubic.py** — exhaustive generator for connected cubic graphs.
- **src/signings/theorems.py** — executable transcriptions of the proof steps
  (bottleneck formula, the pencil $(M_\ell,K_\ell)$, the $(P,R)$ determinant,
  the $\ell=5$ witness, the explicit three-arm construction).
- **tests/test_project.py** — 18 tests, all passing.
- **experiments/** — exploration scripts, including the two that reproduce the
  reported numbers (`run_all_experiments.py`, `case_b_symbolic.py`).

### Documentation

- **README.md** — problem, main theorem, novelty qualification, repository map,
  reproduction commands, limitations, open questions.
- **docs/reproducibility.md** — versions, commands, expected output, and a
  precise statement of what is proved versus what is machine-checked.
- **LICENSE** — MIT.