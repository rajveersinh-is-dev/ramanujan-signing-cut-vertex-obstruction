# Definitions and notation

Throughout, $G$ is a **finite simple undirected graph**.  $A(G)$ is its
$0/1$ adjacency matrix ($A(G)_{uv}=1$ iff $\{u,v\}\in E(G)$, zero diagonal).
All vector spaces are over $\mathbb R$ unless stated otherwise, and all inner
products and norms are Euclidean.

## 1. Signings

A **signing** of $G$ is a map $\sigma\colon E(G)\to\{-1,+1\}$.  Its **signed
adjacency matrix** is
$$
A_{\sigma}(G)_{uv}=\begin{cases}\sigma(\{u,v\}),& \{u,v\}\in E(G),\\ 0,&\text{otherwise},\end{cases}
$$
with zero diagonal.  The **vertex sums** of $\sigma$ are
$$
\sigma(v)=\sum_{e\ni v}\sigma(e)\qquad (v\in V(G)).
$$

**Vertex switching.**  For $\delta\in\{-1,+1\}^{V(G)}$ the *switched* signing is
$(\delta\cdot\sigma)(\{u,v\})=\delta(u)\sigma(\{u,v\})\delta(v)$.  With
$D=\operatorname{diag}(\delta)$ one has
$$
A_{\delta\cdot\sigma}(G)=D\,A_\sigma(G)\,D,
$$
so the spectrum of $A_\sigma$ is a switching invariant.  The orbit of a signing
under switching has size $2^{n-1}$ for $n=|V(G)|$ (the map $\delta\mapsto-\delta$
acts trivially, and no other non-trivial $\delta$ does).  Consequently, for
connected $G$ on $n$ vertices with $m$ edges, there are exactly $2^{m-n+1}$
switching classes, and a set of representatives is obtained by switching along
a spanning tree: this fixes the signs of the $n-1$ tree edges and leaves
$m-n+1$ signs free.

**Remark (a common error).**  Switching does *not* preserve the vertex sums:
$\sigma'(v)=\delta(v)\sum_{w\sim v}\sigma(\{u,w\})\delta(w)$, which in general
differs from $\pm\sigma(v)$.  An earlier draft of this project asserted the
opposite; the mistake was caught by `tests/test_project.py` and is recorded in
`research/proof_attempts.md`.

## 2. Spectral quantities

For $G$ regular of degree $d$:

* $\|A_\sigma\|_{\mathrm{op}}$ is the Euclidean operator norm; since
  $A_\sigma$ is symmetric, $\|A_\sigma\|_{\mathrm{op}}=\max_i|\lambda_i|$ where
  $\lambda_1\ge\lambda_2\ge\cdots\ge\lambda_n$ are its eigenvalues.
* $\rho_d:=2\sqrt{d-1}$ is the **Ramanujan bound** for degree $d$.
* $G$ is **Ramanujan** if every eigenvalue of $A(G)$ other than the trivial
  ones ($\lambda_1=d$, and $\lambda_n=-d$ when $G$ is connected and bipartite)
  lies in $[-\rho_d,\rho_d]$.
* For connected regular $G$ we write $\lambda_2(A(G))$ for the second largest
  eigenvalue of $A(G)$.  Since $G$ is connected, $\mathbf 1$ spans the top
  eigenspace of $A(G)$ and
  $$
  \lambda_2(A(G))=\max\Big\{\tfrac{x^{\mathsf T}A(G)x}{x^{\mathsf T}x}: 0\ne x,\ x\perp \mathbf 1\Big\}.
  \tag{D}
  $$

## 3. The signing norm $\beta$

$$
\beta(G):=\min_{\sigma\colon E(G)\to\{\pm1\}}\ \|A_\sigma(G)\|_{\mathrm{op}}.
\tag{1}
$$
$\beta(G)$ is the best possible spectral radius achievable by signing the edges.

**Note.**  $\max_v|\sigma(v)|$ does *not* upper-bound $\|A_\sigma(G)\|_{\mathrm{op}}$.
For $K_4$ signed with $\sigma(\{1,2\})=\sigma(\{3,4\})=-1$ and all other edges
$+1$, every vertex sum is $\pm1$ while the eigenvalues of $A_\sigma$ are
$3,-3,1,1$, so $\|A_\sigma\|=3$.  The true optimum is
$\beta(K_4)=\sqrt5$ (verified in `tests/test_project.py`).

## 4. The cut-vertex set-up

Let $G$ be a **cubic** graph (3-regular) and $z\in V(G)$ a cut vertex.  Write
$$
G-z=C_1\;\dot\cup\;C_2\;\dot\cup\;C_3,\qquad n_i:=|C_i|,\quad n_1\le n_2\le n_3.
$$

**Lemma 0 (structure at a cut vertex).**  $G$ is connected and cubic, $z$ is a
cut vertex, then each $C_i$ is connected, $|C_i|$ is odd and $\ge 5$, and each
$C_i$ is joined to $z$ by exactly one edge.

*Proof.*  Each component of $G-z$ must contain a neighbour of $z$ (otherwise
$G$ would not be connected), and $z$ has exactly three neighbours, so $t=3$
components and each is joined to $z$ by exactly one edge.  Counting degrees,
$3n_i=2|E(C_i)|+1$, so $n_i$ is odd and $|E(C_i)|=(3n_i-1)/2$.  Since
$|E(C_i)|\le\binom{n_i}{2}$, we get $n_i^2-4n_i+1\ge0$, hence $n_i\ge4$, and
being odd, $n_i\ge5$. $\square$

Denote by $t_i$ the unique vertex of $C_i$ adjacent to $z$.  Two useful families
of vectors (all orthogonal to $\mathbf 1$):

$$
p_i:=\mathbf 1_{C_i}-n_i e_z,\qquad
q_i:=e_{t_i}-\tfrac1{n_i}\mathbf 1_{C_i}.
\tag{2}
$$
Indeed $\langle p_i,\mathbf 1\rangle = n_i-n_i=0$ and
$\langle q_i,\mathbf 1\rangle = 1-\frac1{n_i}n_i=0$.