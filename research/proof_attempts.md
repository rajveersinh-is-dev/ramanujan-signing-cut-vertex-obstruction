# Proof attempts, dead ends, and corrected errors

This file records what was tried and what failed, including claims that were
made during the investigation and later **withdrawn**.  It is deliberately
kept in the repository.

---

## A. Approaches that failed to reach the conjecture

### A1. Perfect matchings (dead end)

**Idea.** For a cubic graph $G$, a signing with all vertex sums $\sigma(v)=\pm1$
is exactly the signing that makes the $-1$ edges a perfect matching.  If such a
signing had $\|A_\sigma\|\le2\sqrt2$ we would be done for all cubic graphs with
a perfect matching.

**Outcome: dead end, and the motivating lemma is FALSE.**  Vertex sums being
small does *not* control the operator norm.  Explicit counterexample computed in
this project: $K_4$ with $\sigma(\{1,2\})=\sigma(\{3,4\})=-1$ has all vertex
sums $\pm1$ but $\|A_\sigma\|=3$ (eigenvalues $3,-3,1,1$).  The true optimum
$\beta(K_4)=\sqrt5$ requires signing with *larger* vertex sums.  The perfect
matching route was abandoned.

### A2. Random signings / second moment (not pursued to completion)

**Idea.** $A_\sigma^2$ has diagonal $d$ and off-diagonal entries
$\sum_{w\in N(u)\cap N(v)}\sigma(uw)\sigma(vw)$, so $|A_\sigma^2|\le A^2$
entrywise and, by the sandwiching principle, $\lambda_{\max}(A_\sigma^2)\le
\max_v\sigma(v)^2+d$ … which reduces to the trivial $\|A_\sigma\|\le d$.  A
Friedman-style second-moment argument would give
$\|A_\sigma^2\|\lesssim\lambda_2(A)^2+d$ with fluctuations $O(d\sqrt{\log n})$,
which beats $\rho_d^2=4(d-1)$ only when $\lambda_2^2+d+(\text{fluctuation})
<4(d-1)$, i.e. only for graphs that are *strictly* Ramanujan with $n$ bounded
in terms of the spectral gap.

**Outcome: not pursued.**  It would give a quantitative statement for small $n$
or near-Ramanujan graphs, not a theorem about all Ramanujan graphs, and would
require a second-moment bound that is a research project in itself.  Recorded as
a possible future direction.

### A3. Exhaustive enumeration of cubic graphs (dropped as a *primary* method)

**Idea.** Enumerate all connected cubic graphs on $n$ vertices and compute
$\beta(G)$ exactly by exhaustive search over switching classes ($2^{n/2+1}$ of
them).

**Outcome: partially implemented, then demoted.**  The generator
`src/signings/cubic.py` is correct (validated against A002851 for $n=4,6,8$)
but the isomorphism filter had to be replaced because `igraph`'s
`canonical_permutation()` is *not* a true canonical labelling (verified: two
isomorphic labelings of $K_{3,3}$ receive different canonical permutations; see
`experiments/_gen_check.py` history).  This would have to be replaced by a
refinement/individualisation canonical form.  Even then, this would be a
*census*, not a theorem, and would dominate the project's runtime.  Demoted to a
possible future computation.

---

## B. Errors made and corrected during the investigation

These are recorded because a reader checking the repository should know which
statements were wrong at some point and how the error was detected.

### B1. "Vertex sums are invariant under switching" — **FALSE, withdrawn**

**Claim made.**  Switching $\sigma\mapsto\delta\cdot\sigma$ "preserves the
spectrum and leaves every vertex sum $\sigma(v)$ unchanged".

**Why it is false.**  $A_{\delta\cdot\sigma}=DA_\sigma D$ so the *spectrum* is
preserved, but $A_\sigma\mathbf 1=\sigma(v)$ gives
$A_{\delta\cdot\sigma}\mathbf 1=DA_\sigma(\delta\mathbf 1)$, whose entries are
$\delta(v)\sum_{w\sim v}\sigma(\{v,w\})\delta(w)$ — this is $\delta(v)\sigma(v)$
only in special cases.  Concrete counterexample: on $C_4$ with
$\sigma\equiv+1$ (all vertex sums $+1$), switching by
$\delta=(1,1,-1,-1)$ gives all vertex sums $0$.

**How detected.**  `tests/test_project.py::test_switching_invariance` failed
even after the spectrum assertion passed.

**Correction.**  The claim was removed from `mathematics/definitions.md` and
from the module docstring of `src/signings/core.py`; only spectral invariance is
asserted now, and the fact is documented as a common error.

### B2. Three hand-derived bilinear forms for the cut-vertex pencil — **FALSE, withdrawn**

During the derivation of Theorem 3 the bilinear forms were computed by hand
three times, and each time produced a pencil giving a Rayleigh quotient larger
than $3$ — impossible for a cubic graph, since $A(G)$ has all eigenvalues in
$[-3,3]$.  Three separate errors were made:

1. The $z$-edge contribution was counted once instead of twice
   ($\sum_i [z_1(c_2{+}d_2)+(c_1{+}d_1)z_2]$, not $2z_1\sum + 2z_2\sum$).
2. The "tip" value was encoded as $c_i=0,d_i=1$ for $q_i$ instead of
   $c_i=-1/n_i,d_i=1$.
3. The arm-internal adjacency form was written as
   $2(n-1)c_1c_2+(c_1+d_1)(c_2+d_2)+(nc_1+d_1)(c_2+d_2)$, which is wrong.

**How detected.**  The computed top eigenvalue exceeded $3$, violating the
spectral bound for a cubic graph; this was caught by the *cross-check against
explicitly constructed coordinate vectors*, which is recorded in
`tests/test_project.py::test_residual_pencil_matches_explicit_vectors`.

**Correction.**  The correct formula
$$x^{\mathsf T}A(G)y\big|_{C_i}=3(n_i-1)c_1c_2+2d_1c_2+2c_1(c_2+d_2)$$
is proved in `mathematics/proofs.md` (item 2 of the derivation), and
`experiments/case_b_symbolic.py` now derives the pencil from the three
identities and is checked against explicit vectors for $\ell=5,7,11,21,101,1001$.

### B3. `igraph.canonical_permutation` used as an isomorphism test — **FALSE, withdrawn**

Verified experimentally: two isomorphic labelings of $K_{3,3}$ get canonical
permutations `[5,4,2,1,3,0]` and `[5,4,1,3,2,0]`.  The function returns *a*
permutation into canonical position, not a canonical form.  The cubic generator
was not shipped with this bug.

### B4. "$\max_v|\sigma(v)|$ bounds $\|A_\sigma\|$" — **FALSE, withdrawn**

A common claim in the literature's heuristic discussion.  Refuted by the $K_4$
example above ($K_4$: all $|\sigma(v)|=1$, $\|A_\sigma\|=3$).  See A1.

### B5. Two false versions of the 2-lift proposition (Result 4) — **WITHDRAWN, corrected**

I first claimed:

> **P1 (false).** If $G$ is a connected Ramanujan 2-lift of $G_0$, then
> $\beta(G_0)\le\rho_d$.

and then the "refinement"

> **P2 (false).** The exceptional value $\varepsilon d\in\mathrm{Spec}(A_{\sigma_0}(G_0))$
> occurs **iff** $G_0$ is non-bipartite.

**Why P1 is false.**  Take $G_0=K_4$ and $\sigma_0\equiv-1$.  Then
$\mathrm{Spec}(A_{\sigma_0}(K_4))=\{-3,1,1,1\}$, so $\|A_{\sigma_0}(K_4)\|=3>2\sqrt2$.
The associated lift is the bipartite double cover of $K_4$: it is connected,
cubic, and has spectrum $\{-3,-1,-1,-1,1,1,1,3\}$, hence it **is** Ramanujan.
Exhaustive enumeration over the five base graphs used found **552** such
counterexample lifts.

**Why P2 is false.**  For $G_0=K_4$ only $8$ of the $64$ signings put $\pm3$ in
$\mathrm{Spec}(A_{\sigma_0}(K_4))$, yet $K_4$ is not bipartite; over the bases
tested there were **32240** violations.

**Underlying error.**  Both used the claim "a 2-lift of $G_0$ is bipartite iff
$G_0$ is".  That is false: the bipartite double cover of a *non-bipartite* cubic
graph is connected, so a lift can be bipartite while its base is not.  The
correct and tautological statement is $-d\in\mathrm{Spec}(A_{\sigma_0}(G_0))$
iff $G$ is bipartite.

**Correction (now Proposition 4 in `mathematics/proofs.md`).**  $G_0$ is
connected; $d$ never occurs in $\mathrm{Spec}(A_{\sigma_0}(G_0))$;
$\mathrm{Spec}(A_{\sigma_0}(G_0))\subseteq[-\rho_d,\rho_d]\cup\{-d\}$; and $-d$
cannot occur when $G_0$ is bipartite.  Hence $\beta(G_0)\le\rho_d$ whenever
$G_0$ is bipartite.  Verified exhaustively: `experiments/verify_lift_proposition.py`
enumerates all $37952$ signings of $K_4$, $K_{3,3}$, the triangular prism, the
cube and the Petersen graph ($37240$ connected lifts, all Ramanujan) with
**0** failures on each of the four checks.

**Lesson.**  The failing step was a *bipartiteness* claim about lifts, i.e. a
statement about a different invariant than the one under study; it would not
have been caught by checking the eigenvalues alone.  Writing the exhaustive
enumerator *before* trusting the proposition is what exposed it.

---

## C. Attempts that succeeded

* **The bottleneck vector for a pair of arms.**  Successful and sharp; became
  Theorem 1.
* **Using the third arm's indicator with a nonzero value at $z$.**  Allowing
  $x_z\ne0$ (vector $\mathbf 1_{C_i}-n_ie_z$ instead of $\frac1{n_i}\mathbf
  1_{C_i}$) was the key improvement over the naive bottleneck vector: it raised
  the bound for arm orders $(5,5,\ell)$ from $2.800$ to $2.851$ at $\ell=5$.
* **Adding a mean-zero "tip-dip" direction.**  The vectors
  $q_i=e_{t_i}-\frac1{n_i}\mathbf 1_{C_i}$ supply the extra $0.051$ needed in
  the residual case; became Theorem 3.
* **Reducing to a single arm-order family.**  Corollary 1 shows all
  configurations with $n_2\ge7$ are excluded by the plain bottleneck bound, so
  only the one-parameter family $(5,5,\ell)$ needed the richer test space.  This
  made the symbolic computation tractable.