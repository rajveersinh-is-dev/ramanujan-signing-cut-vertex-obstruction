# Literature review and known results

All references below were verified to exist during the present investigation
(Oct 2026).  Where a statement about "current status" is made, the source and
date are given.

## 1. The signing problem

**Bilu and Linial**, *Lifts, discrepancy and nearly optimal spectral gap*,
Combinatorica **26** (2006) 495–519.  Asked whether every finite simple
$d$-regular graph admits a signing $\sigma$ with
$\|A_\sigma(G)\|\le 2\sqrt{d-1}$.  Xu's paper cites this as Conjecture 3.1.

**Bilu and Linial**, *Ramanujan signing of regular graphs*, Combin. Probab.
Comput. **13** (2004) 911–912.  An earlier, weaker question.

**Marcus, Spielman and Srivastava**, *Interlacing families I: bipartite Ramanujan
graphs of all degrees*, Ann. of Math. (2) **182** (2015) 307–325, arXiv:1304.4132.
Theorem 5.3 of that paper states that every $d$-regular graph ($d\ge2$) admits a
signing whose signed adjacency eigenvalues are **all at most** $2\sqrt{d-1}$ —
a *one-sided* statement.  As Xu (2026, §1.3) points out, this settles the
Bilu–Linial conjecture only for bipartite graphs, where the signed spectrum is
symmetric about $0$.  This distinction (upper bound vs. operator norm) is the
crux of the current problem, and it is why MSS does **not** contradict Xu.

**Ravichandran and Srivastava**, *Asymptotically optimal multi-paving*,
IMRN 2021, 10908–10940.  Bound $2\sqrt{2(d-1)}$ on the norm.

**Hall, Puder and Sawin**, *Ramanujan coverings of graphs*, Adv. Math. **323**
(2018) 367–397.  Every connected $d$-regular graph ($d\ge3$) has an $\ell$-lift
with new eigenvalues $\le2\sqrt{d-1}$ for $\ell\ge2$.

**Zhiqiang Xu**, *A 3-regular counterexample to the Bilu--Linial signing
conjecture*, arXiv:2609.15591v4 (submitted 14 Sep 2026, revised 29 Sep 2026).
The abstract states: "We construct a finite connected simple cubic graph $F$
such that every signing of its edges yields a signed adjacency matrix with an
eigenvalue outside $[-2\sqrt2,2\sqrt2]$. This disproves the Bilu--Linial signing
conjecture for general regular graphs. The graph $F$ is not Ramanujan, and the
conjecture restricted to Ramanujan base graphs remains open."

  *Verification:* fetched from `https://arxiv.org/abs/2609.15591` and the full
  text `https://arxiv.org/html/2609.15591v4` on 12–13 Oct 2026.

**Xu, Section 8 (Conclusion), verbatim paraphrase.** Two questions remain open:
(a) "does every finite connected simple regular *Ramanujan* graph admit a
signing whose entire spectrum lies in the Ramanujan interval?" — this is
Conjecture R of this project; (b) "does every finite connected simple regular
graph admit an $\ell$-lift, for some integer $\ell\ge2$, whose new eigenvalues
all lie in the Ramanujan interval?"

**Related recent work** (listed by Xu, §1.3, and by arXiv listing):

| Reference | Statement |
|---|---|
| Huang, *A $\sqrt2$-approximation to the Bilu--Linial conjecture*, arXiv:2609.17551 (2026) | graph-theoretic reduction to the bipartite signing theorem |
| Lin and Zhou, arXiv:2609.15715 (2026) | $\beta(G)<\tfrac{3+\sqrt5}{2}\sqrt{d-1}$ for all $d$-regular $G$ |
| Xu, Xu and Zhu, *Improved bounds in Weaver's $\mathrm{KS}_r$ conjecture*, J. Funct. Anal. **285** (2023) 109978 | higher-rank interlacing families |

**Status: the Bilu--Linial conjecture is FALSE (disproved 2026-09).
Conjecture R (Ramanujan signing) is OPEN as of 2026-10-13.**

**Priority check (re-run 2026-10-13).**  The arXiv API query `all:"Bilu-Linial"`
returns 14 entries, all inspected.  The most recent are

| arXiv | Date | Title | Bears on Conjecture R? |
|---|---|---|---|
| 2609.23855 | 2026-09-20 | *Kadison--Singer partitions and Bilu--Linial graph signings in polynomial time* | no --- algorithmic |
| 2609.15715 | 2026-09-14 | *Improved bounds ... via spectral recovery* (Lin--Zhou) | no --- weaker upper bound |
| 2609.15591 | 2026-09-14, v4 09-29 | *A 3-regular counterexample ...* (Xu) | **states Conjecture R as open** |
| 2609.17551 | 2026-07-16 | *A $\sqrt2$-approximation to the Bilu--Linial conjecture* (Huang) | no --- weaker upper bound |
| 2606.28797 | 2026-06-27 | *An improved upper bound ... via interlacing families* | no --- weaker upper bound |

None of these addresses the Ramanujan-signing question; the frontier is where
Xu's §8 leaves it.  Several of them still phrase the Bilu--Linial conjecture as
open, which reflects preprint lag rather than disagreement with Xu.

## 2. What is standard background used here

* **Bottleneck / "two-lobe" Rayleigh bound for regular graphs.**  For a
  $d$-regular graph and a vertex subset $S$ with $c=e(S,S^c)$ the test vector
  $\alpha\mathbf 1_S+\beta\mathbf 1_{S^c}$ orthogonal to $\mathbf 1$ gives
  $\lambda_2\ge d-c\bigl(\tfrac1{|S|}+\tfrac1{n-|S|}\bigr)$.  This is a
  standard consequence of the variational characterisation of $\lambda_2$; we
  use a sharper *two-arm* version (Theorem 1) which balances two components of
  $G-z$ rather than $S$ against its whole complement.
* **Structure of a cubic graph at a cut vertex.**  Classical: $G-z$ has exactly
  three components, each of odd order $\ge5$, each attached to $z$ by exactly
  one edge.  Proved from scratch here as Lemma 0.
* **Petersen's theorem / perfect matchings.**  Considered and *not* used: for a
  cubic graph, signings with $|\sigma(v)|\le1$ are exactly the perfect-matching
  signings, but as shown in `mathematics/definitions.md` a small vertex-sum does
  **not** imply a small operator norm, so this route gives nothing.
* **A002851** (number of connected cubic graphs on $2n$ vertices: 1, 2, 5, 19,
  85, 509, 4060, …) is used only to validate the generator in
  `src/signings/cubic.py`.

## 3. Known status of things this project could have attacked but did not

* **Smallest cubic graph with $\beta>2\sqrt2$.**  Xu's example has an enormous
  number of vertices.  Whether a small one exists is not addressed in the
  literature located here and is not settled by this project.  Given
  $2^{n/2+1}$ switching classes per graph and 4060 connected cubic graphs on 16
  vertices, exhaustive search to $n=16$ is computationally trivial, so this is
  listed as a feasible future computation, not a mathematical obstacle.
* **Conjecture R for $d\ge4$.**  Nothing found.

## 4. Searches performed

* arXiv API queries for `"Bilu-Linial"`, `ti:"Zero-sum Ramsey"`,
  `all:"cross number" AND all:"zero-sumfree"`, `abs:"remains open" AND cat:math.CO`,
  `abs:"is open" AND cat:math.NT` (Oct 2026).
* Crossref bibliographic queries for the zero-sum Ramsey and little-cross-number
  candidates.
* Semantic Scholar graph search (returned empty/rate-limited for this topic).
* Direct fetch of `arxiv.org/abs/2609.15591` and its HTML full text.
* DuckDuckGo HTML search for "Ringel-Jungerman" and related terms (no results;
  that candidate was dropped).