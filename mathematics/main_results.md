# Main results

All statements below are **PROVED HERE** unless otherwise labelled.  Complete
proofs are in `proofs.md`; the same proofs are mirrored in
`paper/main.tex`.

---

## Result 1 (PROVED HERE): the two-arm bottleneck bound

**Theorem 1.** Let $G$ be a connected cubic graph and $z$ a cut vertex of $G$.
Write $G-z=C_1\dot\cup C_2\dot\cup C_3$ with $n_i=|C_i|$ and $n_1\le n_2\le n_3$.
Then
$$
\lambda_2(A(G))\ \ge\ \max_{i<j}\ \Big(3-\frac{n_i^2+n_j^2}{n_in_j(n_i+n_j)}\Big).
\tag{3}
$$
In particular
$$
\lambda_2(A(G))\ \ge\ \max_{i<j}\Big(3-\frac{1}{\max(n_i,n_j)}\Big).
\tag{4}
$$

**Comment.** The bound (3) is *sharp in order*: for arms of orders $n$ and $cn$
with $c$ fixed, (3) tends to $3-c^{-1}n^{-1}\to 3$, which is the true limit of
$\lambda_2$ for three-arm graphs whose arms all grow (see
`experiments/cut_vertex_spectrum.py`).  The estimate itself — the Rayleigh
quotient of a two-component indicator vector — is classical ("bottleneck"
or "Cheeger-type" bound); what is proved here is the statement in this exact
form together with its use in Result 2.

---

## Result 2 (PROVED HERE): Ramanujan cubic graphs have no cut vertex

**Corollary 1.** If $G$ is connected and cubic and $z$ is a cut vertex of $G$
such that at least two of the three components of $G-z$ have order $\ge 7$, then
$$
\lambda_2(A(G))\ \ge\ \frac{20}{7}\ >\ 2\sqrt2 .
\tag{5}
$$

**Theorem 2.** *Every connected cubic graph with a cut vertex is not Ramanujan*:
$$
\lambda_2(A(G))>2\sqrt2 .
\tag{6}
$$

**Corollary 2.** Every connected Ramanujan cubic graph is 2-connected; in
particular $G-z$ being empty of cut vertices forces any cut vertex of $G$ to
split $G$ into three arms of orders $\ge5$ with at least two of them equal to $5$
(the only configuration not excluded by (5) alone).

**Why this matters.** The conjecture under study (Conjecture R below) asserts
that Ramanujan graphs admit Ramanujan signings.  Theorem 2 shows that for
$d=3$ any counterexample must live among graphs that are at least 2-connected —
in particular Xu's (2026) counterexample, whose mechanism is a cut vertex with
three huge arms, cannot be reused against Conjecture R.  It also shows that the
one-sided bound of Marcus–Spielman–Srivastava is the only known general
statement, and that the two-sided question is genuinely harder.

---

## Result 3 (PROVED HERE): the residual four-dimensional test space

**Theorem 3.** Suppose $n_1=n_2=5$ and $n_3=\ell\ge5$ (odd).  Put
$$
W=\operatorname{span}\{\,p_1+p_2,\ q_1+q_2,\ p_3,\ q_3\,\}\subseteq\mathbb R^{V(G)}\cap\mathbf 1^\perp .
$$
Then
$$
\lambda_2(A(G))\ \ge\ \lambda^*(\ell):=\lambda_{\max}\big(M_\ell,K_\ell\big)\ >\ 2\sqrt2
\quad\text{for every odd }\ell\ge5,
\tag{7}
$$
where $M_\ell= p^{\mathsf T}A(G)p$ and $K_\ell=p^{\mathsf T}p$ for the
$4\times 4$ matrix $p$ whose columns are the displayed generators, namely

$$
M_\ell=\begin{pmatrix}
-12 & -\tfrac{88}{5} & -2\ell-10 & -10+\tfrac{10}{\ell}\\[2pt]
-\tfrac{88}{5} & -\tfrac{12}{25} & -\tfrac{8\ell}{5} & 0\\[2pt]
-2\ell-10 & -\tfrac{8\ell}{5} & \ell-1 & -\ell+\tfrac1\ell\\[2pt]
-10+\tfrac{10}{\ell} & 0 & -\ell+\tfrac1\ell & -\tfrac1\ell-\tfrac1{\ell^2}
\end{pmatrix},\quad
K_\ell=\begin{pmatrix}
110 & 0 & 10\ell & 0\\
0 & \tfrac85 & 0 & 0\\
10\ell & 0 & \ell^2+\ell & 0\\
0&0&0& 1-\tfrac1\ell
\end{pmatrix}.
$$

Numerically $\lambda^*(5)=2.850781$, and $\lambda^*(\ell)\nearrow2.95614$ as
$\ell\to\infty$; the minimum over odd $\ell\ge5$ occurs at $\ell=5$.

The pencil depends only on the *arm orders*, not on the internal structure of
the arms, so (7) is a genuine uniform statement over all cubic graphs with a cut
vertex whose two smallest arms have order $5$.

---

## The open problem (not solved here)

**Conjecture R (Xu, 2026; open).** *Ramanujan signing.* If $G$ is a finite
connected simple $d$-regular Ramanujan graph ($d\ge2$), there is a signing
$\sigma$ with $\|A_\sigma(G)\|_{\mathrm{op}}\le 2\sqrt{d-1}$.

**Conjecture R3 (its $d=3$ case; open).** If $G$ is a connected cubic graph with
$\lambda_2(A(G))\le 2\sqrt2$, then $\beta(G)\le2\sqrt2$.

Neither is proved or disproved in this repository.  See `research/known_results.md`
for the exact status evidence and `research/proof_attempts.md` for the attacks
tried.

---

## Evidence labels used in this project

| Label | Meaning |
|---|---|
| **PROVED HERE** | complete proof in `mathematics/proofs.md`, audited |
| **ESTABLISHED IN PRIOR WORK** | supported by a verified citation |
| **COMPUTATIONALLY VERIFIED** | checked on an explicitly described finite domain |
| **CONJECTURED** | evidence but no proof |
| **UNRESOLVED** | an essential step remains open |
| **HEURISTIC** | directional only |