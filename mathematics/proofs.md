# Proofs

Notation is as in `definitions.md`.  Every claim is proved; no step is left as
"clearly" or "obvious".

---

## Lemma 0 (structure at a cut vertex)

*Let $G$ be connected and cubic and $z$ a cut vertex.  Then $G-z$ has exactly
three components $C_1,C_2,C_3$, each of odd order $\ge5$, each joined to $z$ by
exactly one edge; $t:=3$.*

**Proof.** Let $C_1,\dots,C_t$ be the components of $G-z$ and
$e_i:=|\{\text{edges between }C_i\text{ and }z\}|$.  Every $e_i\ge1$, otherwise
$C_i\cup\{z\}$ would be a union of components of $G$ separate from the rest,
contradicting connectedness of $G$ (since $z$ is incident to at least one edge
into each other $C_j$, $j\ne i$).  Also $\sum_i e_i=\deg(z)=3$, hence $t=3$ and
$e_i=1$ for every $i$.

Counting degrees inside $C_i$: each of its $n_i$ vertices has degree $3$, and
exactly one of those edges leaves $C_i$, so $3n_i=2|E(C_i)|+1$.  Thus $n_i$ is
odd and $|E(C_i)|=(3n_i-1)/2$.  As $C_i$ is a simple graph,
$\tfrac{3n_i-1}{2}\le\tfrac{n_i(n_i-1)}{2}$, i.e. $n_i^2-4n_i+1\ge0$.  The
roots are $2\pm\sqrt3$, so $n_i\le 2-\sqrt3<1$ or $n_i\ge 2+\sqrt3>3.7$; hence
$n_i\ge4$, and $n_i$ odd gives $n_i\ge5$. $\square$

---

## Theorem 1 (two-arm bottleneck bound)

*With the hypotheses of Lemma 0 and $n_1\le n_2\le n_3$,
$\lambda_2(A(G))\ge\max_{i<j}\big(3-\frac{n_i^2+n_j^2}{n_in_j(n_i+n_j)}\big)$.*

**Proof.** Fix two distinct components $C_i,C_j$ and write $a=n_i$, $b=n_j$.
Consider
$$
x=\tfrac1a\,\mathbf 1_{C_i}-\tfrac1b\,\mathbf 1_{C_j}.
$$
Then $\langle x,\mathbf 1\rangle=\tfrac a a-\tfrac b b=0$, so by (D) it suffices
to compute $R(x)=x^{\mathsf T}A(G)x/\|x\|^2$.

*Norm.*  $\|x\|^2=\tfrac1{a^2}\cdot a+\tfrac1{b^2}\cdot b=\tfrac1a+\tfrac1b$,
because $C_i$ and $C_j$ are disjoint and $x$ vanishes outside $C_i\dot\cup C_j$.

*Quadratic form.*  Since $A$ is symmetric with zero diagonal,
$2x^{\mathsf T}Ax=2\sum_{\{p,q\}\in E}x_p x_q$.  There are no edges between
$C_i$ and $C_j$; edges leaving $C_i$ go only to $z$ where $x_z=0$; likewise for
$C_j$.  Hence
$$
x^{\mathsf T}Ax=\tfrac1{a^2}\,2|E(C_i)|+\tfrac1{b^2}\,2|E(C_j)|
=\tfrac{3a-1}{a^2}+\tfrac{3b-1}{b^2},
$$
using $2|E(C_k)|=3n_k-1$ from Lemma 0.

Putting the numerator over the common denominator $a^2b^2$ and the
denominator over $ab$,
$$
R(x)=\frac{b^{2}(3a-1)+a^{2}(3b-1)}{ab(a+b)}
=\frac{3ab^{2}-b^{2}+3a^{2}b-a^{2}}{ab(a+b)}
=3-\frac{a^{2}+b^{2}}{ab(a+b)} ,
$$
the last step being the identity
$3ab^{2}+3a^{2}b-a^{2}-b^{2}=3ab(a+b)-(a^{2}+b^{2})$, which holds term by term.

Since $R(x)$ is one of the quotients maximised in (D),
$\lambda_2(A(G))\ge R(x)$ for this pair, and maximising over the three pairs
gives (3). $\square$

---

---

## Corollary 1

*If at least two of $n_1,n_2,n_3$ are $\ge 7$, then $\lambda_2(A(G))\ge 20/7>2\sqrt2$.*

**Proof.** Write $a\le b$ for two arm orders with $a\ge7$.  We need an upper
bound on $\frac{a^2+b^2}{ab(a+b)}$.  Put $t=a/b\in(0,1]$.  Then
$$
\frac{a^2+b^2}{ab(a+b)}=\frac{t^2+1}{tb(t+1)}=\frac1a\cdot\psi(t),
\qquad \psi(t):=\frac{t+1/t}{1+t}.
$$
For $0<t\le1$,
$$
\psi(t)\le1 \iff t+\tfrac1t\le 1+t \iff \tfrac1t\le 1 ,
$$
which holds.  Hence $\frac{a^2+b^2}{ab(a+b)}\le\frac1a\le\frac17$, and by (3)
$$
\lambda_2(A(G))\ge 3-\tfrac17=\tfrac{20}{7}.
$$
Finally $\frac{20}{7}>2\sqrt2$, because $\frac{400}{49}>8\iff400>392$. $\square$

---

## Theorem 2 (cubic graphs with a cut vertex are never Ramanujan)

*$\lambda_2(A(G))>2\sqrt2$.*

**Proof.** By Lemma 0 and the ordering $n_1\le n_2\le n_3$, two cases.

**Case 1: $n_2\ge7$.**  Apply Corollary 1 with the pair $(n_2,n_3)$:
$\lambda_2(A(G))\ge\frac{20}{7}>2\sqrt2$.

**Case 2: $n_2=5$.**  Then $n_1=5$ and $n_3=\ell$ is odd with $\ell\ge5$.  We
use the four-dimensional subspace
$$
W=\operatorname{span}\{P,Q,R,S\},\qquad
P=p_1+p_2,\;Q=q_1+q_2,\;R=p_3,\;S=q_3
$$
with $p_i,q_i$ as in `definitions.md` (2).  Every element of $W$ is orthogonal
to $\mathbf 1$, so $\lambda_2(A(G))\ge\lambda_{\max}(M_\ell,K_\ell)$ with
$M_\ell=p^{\mathsf T}A(G)p$, $K_\ell=p^{\mathsf T}p$.

**Computation of the pencil.**  We record the three combinatorial identities
used, and prove them.  For a vector $x$ supported on $\{z\}\cup C_1\cup C_2\cup
C_3$ let $(z_x)$ be its value at $z$, and for each arm $i$ write $c_i^x$ for
its value on $C_i\setminus\{t_i\}$ and $d_i^x$ for the *excess* at the tip, so
that $x(t_i)=c_i^x+d_i^x$.

1. **Gram.**  On one arm,
   $\sum_{u\in C_i}x(u)y(u)=(n_i-1)c^xc^y+(c^x+d^x)(c^y+d^y)$, because $C_i$
   has $n_i-1$ non-tip vertices and one tip.  Summing over the three arms and
   adding $z_xz_y$ gives (Gr).

2. **Arm-internal adjacency form.**  Since every vertex of $C_i$ except $t_i$
   has degree $3$ inside $C_i$ and $t_i$ has degree $2$ there,
   $$
   x^{\mathsf T}A(G)y\big|_{C_i}
   =\sum_{u\in C_i}x(u)\!\!\sum_{\substack{v\in C_i\\ v\sim u}}\!\!y(v).
   $$
   Group the inner sum by its second index $v\ne t_i$: each such $v$ is
   adjacent to $d(v)=3$ vertices of $C_i$, of which at most one is $t_i$; since
   $t_i$ is a single vertex, for $v\ne t_i$ the contribution to
   $\sum_u x(u)[\ldots]$ is $c^x\,(3\,y(v))$ unless $t_i\sim v$, in which case
   the value $x(t_i)$ replaces $c^x$ once, contributing $d^x y(v)$.  Summing:
   $c^x\cdot 3\sum_{v\ne t_i}y(v)+\big(\sum_{v\sim t_i,\,v\ne t_i}y(v)\big)d^x
   =3(n_i-1)c^xc^y+2d^xc^y$, where the factor $2$ is the degree of $t_i$
   inside $C_i$.  Finally the pair $\{t_i,u\}$ is counted a second time (once
   from $u$, once from $t_i$), adding $y(t_i)\sum_{u\sim t_i,\,u\ne t_i}x(u)
   =2(c^y+d^y)c^x$.  Hence
   $$
   x^{\mathsf T}A(G)y\big|_{C_i}=3(n_i-1)c^xc^y+2d^xc^y+2c^x(c^y+d^y).
   \tag{*}
   $$

3. **The three edges $zt_i$.**  Each edge $\{z,t_i\}$ contributes
   $z_x\,y(t_i)+x(t_i)\,z_y$ to $x^{\mathsf T}Ay$.

Substituting the generators
$$
P:\ (z_P=-10;\ \ (1,1,0),(2,1,0)),\quad
Q:\ (z_Q=0;\ \ (1,-\tfrac15,1),(2,-\tfrac15,1)),
$$
$$
R:\ (z_R=-\ell;\ (3,1,0)),\quad
S:\ (z_S=0;\ (3,-\tfrac1\ell,1))
$$
into (Gr), (*) and the $z$-edge rule and collecting gives exactly $M_\ell$ and
$K_\ell$ displayed in `main_results.md` Theorem 3.  This computation is
machine-checked: `experiments/case_b_symbolic.py` derives $M_\ell,K_\ell$ from
the three identities and `tests/test_project.py::test_residual_pencil_matches_explicit_vectors`
compares them with matrices built from explicit coordinate vectors on concrete
three-arm graphs, for $\ell=5,7,9,11,21,101,1001$ (agreement to $10^{-9}$).

**Sub-case 2a: $\ell\ge7$.**  Let $D_\ell=M_\ell-2\sqrt2\,K_\ell$.  The
principal minor of $D_\ell$ on indices $\{0,2\}$ (the pair $P,R$) equals
$$
\det D_\ell[\{0,2\},\{0,2\}]
=-4(\ell+11)\big((14\sqrt2-19)\ell-5\sqrt2+2\big).
\tag{**}
$$
The factor $14\sqrt2-19>0$ (since $392>361$).  The remaining factor is positive
as soon as
$$
\ell>\frac{5\sqrt2-2}{14\sqrt2-19}=\frac{67\sqrt2+102}{31}\approx6.347,
$$
so for every integer $\ell\ge7$ it is positive, hence the minor in (**) is
negative.  A symmetric $2\times2$ matrix with negative determinant is
indefinite, so this principal submatrix of $D_\ell$ admits a vector $a\ne0$ with
$a^{\mathsf T}D_\ell a>0$; for such an $a$,
$$
a^{\mathsf T}M_\ell a>2\sqrt2\,a^{\mathsf T}K_\ell a
\ \Longrightarrow\ 
\frac{a^{\mathsf T}M_\ell a}{a^{\mathsf T}K_\ell a}>2\sqrt2 ,
$$
so $\lambda_{\max}(M_\ell,K_\ell)>2\sqrt2$ and $\lambda_2(A(G))>2\sqrt2$.

**Sub-case 2b: $\ell=5$.**  Take the rational vector
$$
a=\Big(-\tfrac35,\ \tfrac15,\ \tfrac65,\ -\tfrac25\Big)^{\mathsf T}\in\mathbb Q^4
$$
(the coordinates refer to $P,Q,R,S$).  Direct substitution into the displayed
$M_5,K_5$ gives
$$
a^{\mathsf T}M_5a-2\sqrt2\,a^{\mathsf T}K_5a
=\frac{19584}{625}-\frac{2748\sqrt2}{125}
=\frac{19584-13740\sqrt2}{625}.
$$
Now $19584^2-2\cdot13740^2=383\,533\,056-377\,515\,200=6\,017\,856>0$, and
$19584>0$, so $19584-13740\sqrt2>0$.  Hence again
$\lambda_2(A(G))\ge\frac{a^{\mathsf T}M_5a}{a^{\mathsf T}K_5a}>2\sqrt2$.

This completes Case 2 and the proof. $\square$

---

## Corollary 2

*Every connected Ramanujan cubic graph is 2-connected.*

**Proof.** By Theorem 2 a connected cubic graph with a cut vertex satisfies
$\lambda_2(A(G))>2\sqrt2=\rho_3$.  Ramanujan-ness for a connected $d$-regular
graph requires *every* non-trivial eigenvalue, in particular the second largest
for a connected graph (whose top eigenspace is $\operatorname{span}\{\mathbf 1\}$),
to have modulus at most $\rho_d$.  That fails. $\square$

---

## Proposition 4 (2-lifts)

*Let $G$ be a finite connected simple $d$-regular Ramanujan graph which is a
2-lift of a $d$-regular graph $G_0$ via the signing $\sigma_0$. Then $G_0$ is
connected, $d\notin\mathrm{Spec}(A_{\sigma_0}(G_0))$,
$\mathrm{Spec}(A_{\sigma_0}(G_0))\subseteq[-\rho_d,\rho_d]\cup\{-d\}$, and $-d$
cannot occur if $G_0$ is bipartite. Hence $\beta(G_0)\le\rho_d$ whenever $G_0$
is bipartite.*

**Setup.** A 2-lift of $G_0$ has $V(G)=V(G_0)\times\{0,1\}$ and is determined by
a signing $\sigma_0$ of $G_0$: an edge $\{u,v\}$ with $\sigma_0(\{u,v\})=+1$
lifts to $\{u^0,v^0\},\{u^1,v^1\}$ and with $\sigma_0(\{u,v\})=-1$ to
$\{u^0,v^1\},\{u^1,v^0\}$.

**The 2-lift spectral identity.** For the uniform signing of the lifted graph,
$$
\mathrm{Spec}\big(A(G)\big)=\mathrm{Spec}\big(A(G_0)\big)\ \sqcup\
\mathrm{Spec}\big(A_{\sigma_0}(G_0)\big).
\tag{L}
$$
(Standard 2-lift identity of Bilu--Linial: decomposing
$\mathbb R^{V(G_0)}\otimes\mathbb R^2$ into the diagonal and anti-diagonal
subspaces, on which $A(G)$ acts as $A(G_0)$ resp. $A_{\sigma_0}(G_0)$. Cited,
not re-proved.)

**Step 1: $G_0$ is connected.** By (L), $\mathrm{Spec}(A(G_0))\subseteq
\mathrm{Spec}(A(G))$ as multisets.  The eigenvalue $d$ occurs exactly once in
$\mathrm{Spec}(A(G))$ because $G$ is connected, and at least once in
$\mathrm{Spec}(A(G_0))$ because $G_0$ is $d$-regular.  Hence it occurs exactly
once in $\mathrm{Spec}(A(G_0))$, so the components of $G_0$ number exactly one:
$G_0$ is connected.

**Step 2: $d$ is in the base part, not the signed part.**  By Step 1 the base
part contains $d$ with multiplicity one, which is the full multiplicity of $d$
in $\mathrm{Spec}(A(G))$.  By (L) the remaining eigenvalues of $A(G)$ — which
are exactly $\mathrm{Spec}(A_{\sigma_0}(G_0))$ — therefore avoid $d$.

**Step 3: the only possible exception is $-d$.**  By Step 1 the eigenvalue $-d$
belongs to $\mathrm{Spec}(A(G))$ if and only if $G$ is bipartite, and then with
multiplicity one.  If $G_0$ is bipartite then $-d\in\mathrm{Spec}(A(G_0))$
(multiplicity one) and, by (L), $-d$ is exhausted by the base part, so it is
absent from the signed part.  If $G_0$ is not bipartite then $-d\notin
\mathrm{Spec}(A(G_0))$, so whenever $-d\in\mathrm{Spec}(A(G))$ (i.e. whenever
$G$ is bipartite) it must lie in the signed part; and conversely if $-d$ is in
the signed part then $G$ is bipartite.

**Conclusion.** Every eigenvalue of $A_{\sigma_0}(G_0)$ is a non-trivial
eigenvalue of $A(G)$, except possibly $-d$; by Ramanujan-ness of $G$ each such
non-trivial eigenvalue lies in $[-\rho_d,\rho_d]$.  By Step 3 the value $d$
never occurs, and if $G_0$ is bipartite then $-d$ does not occur either. $\square$

**Remark.** A lift-compatible signing $\tilde\tau$ of $G$ satisfies
$\mathrm{Spec}(A_{\tilde\tau}(G))=\mathrm{Spec}(A(G_0))\sqcup
\mathrm{Spec}(A_\tau(G_0))$, which contains the eigenvalue $d$.  So a
lift-compatible signing is *never* Ramanujan for $d\ge3$: the 2-lift structure
alone can never produce a Ramanujan signing of $G$, only of its base.

**Correction note — two earlier versions of this proposition were false.**

* *Version 1* claimed $\beta(G_0)\le\rho_d$ for every base graph.  Counterexample:
  $G_0=K_4$, $\sigma_0\equiv-1$.  Then
  $\mathrm{Spec}(A_{\sigma_0}(K_4))=\{-3,1,1,1\}$ (so $\|A_{\sigma_0}(K_4)\|=3>
  2\sqrt2$), while the associated lift — the bipartite double cover of $K_4$ —
  is connected, cubic, and has spectrum $\{-3,-1,-1,-1,1,1,1,3\}$, hence **is**
  Ramanujan.  Exhaustive search found 552 such lifts.
* *Version 2* claimed the exceptional value occurs **iff** $G_0$ is non-bipartite.
  False: for $G_0=K_4$ only $8$ of the $64$ signings give $\pm3$ in
  $\mathrm{Spec}(A_{\sigma_0}(G_0))$, yet $K_4$ is not bipartite;
  $32240$ counterexamples among the bases tested.
* *Underlying error in both*: the claim that a 2-lift of $G_0$ is bipartite iff
  $G_0$ is.  The bipartite double cover of a non-bipartite cubic graph is
  connected, so a lift can be bipartite while its base is not.  The correct
  statement is the tautology $-d\in\mathrm{Spec}(A_{\sigma_0}(G_0))\iff G$ is
  bipartite, which is what Step 3 uses.
* *How detected*: `experiments/verify_lift_proposition.py` enumerates **all**
  signings of $K_4$, $K_{3,3}$, the triangular prism, the cube and the Petersen
  graph ($37952$ signings, $37240$ connected lifts, all of them Ramanujan as
  these graphs are) and checks the claimed conclusions.  Both false versions
  were produced failures.

---

## Proof-audit notes

* The only place where a *finite* computation enters the proof is the identity
  (**) and the two evaluations at $\ell=5$; both are exact and are
  machine-checked against an independent construction.
* The bound (3) is used only in the direction $\lambda_2\ge$ (bound), and the
  proof never assumes the bound is attained.
* Exceptional configurations are handled exhaustively: by Lemma 0 the only
  possible arm orders are odd $\ge5$, and Case 1/Case 2 is a dichotomy on
  $n_2\ge7$ vs $n_2=5$ with $n_1=5$ forced.  There is no third possibility
  because $n_1\le n_2\le n_3$ and $n_1\ge5$.
* Degenerate cases: $z$ must exist as a cut vertex; a graph with no cut vertex
  is outside the scope of Theorems 1–2 and is not needed for Corollary 2.
* **Independent randomised check.**  `experiments/verify_theorem2_random.py`
  builds 4000 random three-arm cubic graphs (random cubic arms, arm orders in
  $5..11$, seed 20261013), computes the true $\lambda_2$, and checks both that
  none is Ramanujan and that the *provable* bound used in the proof never
  exceeds the true $\lambda_2$.  Result: 0 Ramanujan graphs, 0 bound violations,
  smallest margin $\lambda_2-2\sqrt2=0.027345$, attained at arm orders
  $(5,5,5)$ — which is the case the proof treats separately, as predicted.
  This is a consistency check, not part of the proof.
* **Weakest step.**  The step most open to challenge is the arm-adjacency
  formula (*) in the Theorem 3 derivation; it is proved combinatorially in the
  text and independently re-derived symbolically in
  `experiments/case_b_symbolic.py`, then checked against matrices built from
  explicit coordinate vectors on concrete graphs.  If (*) were wrong, the
  resulting bound would (as happened three times during this work) exceed $3$,
  which is impossible for a cubic graph — so the check is sharp.