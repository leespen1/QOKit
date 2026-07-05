# The conditioning correction Var(E[T|c]): derivation and verification

*Written 2026-07-05. Closes (up to a quantified 4% remainder) the open piece of
`research/v2_density_law.md`: Var(E[T|c]), measured ≈60·τ²/m vs the crude 36·τ²/m.
Verified by `research/experiments/015_v2-density-law/conditioning_derivation_check.jl`
(45 instances: ER(0.5), ER(0.25), BA(4), n=12–14; every moment identity below is
asserted as an **exact integer identity**; projections by full 2^n enumeration).*

Notation: spins $s_j=(-1)^{y_j}$, edge products $X_e=s_js_k$,
$S=\sum_{e\in E}X_e = m-2c$ (so conditioning on $c$ ≡ conditioning on $S$),
$T=\sum_i(\sum_{j\sim i}s_j)^2$, codegrees $A_{jk}$, $\tau$ = #triangles,
$q$ = #4-cycles, $P_1=\sum_v\binom{\deg v}{2}$ (edge pairs sharing a vertex).

## Lemma A (quartic split — exact, the centerpiece)

$$\boxed{\,T \;=\; S^2 + m - R_4\,},\qquad
R_4=\sum_{\substack{e\ne f\ \text{ordered}\\ e\cap f=\emptyset}}X_eX_f .$$

**Proof.** $S^2=\sum_e X_e^2+\sum_{e\ne f}X_eX_f = m + \Sigma_{\rm shared} + R_4$,
and for $e=(i,j)$, $f=(i,k)$: $X_eX_f=s_i^2s_js_k=s_js_k$, so
$\Sigma_{\rm shared}=\sum_{j<k}2A_{jk}s_js_k = T-2m$ (density-law Lemma 2). ∎

Consequences. (i) $\mathbb E[T\mid c] = (m-2c)^2 + m - \mathbb E[R_4\mid c]$: an
**exact quadratic-in-$c$ part** — what the crude linear heuristic missed.
(ii) $S^2+m$ is class-constant, so $\operatorname{Var}_{S_v}[T]=
\operatorname{Var}_{S_v}[R_4]$ exactly: **V₂ is the within-class variance of the
disjoint-edge quartic** $R_4$.

## Lemma B (mixed moments — exact, uniform spins)

By monomial cancellation ($\mathbb E[\prod s]{=}1$ iff every vertex appears evenly;
triangles are the only 3-edge even sets, 4-cycles the only 4-edge ones beyond
doubled pairs):
$$\operatorname{Var}(S)=m,\quad \mathbb E[S^3]=6\tau,\quad
\operatorname{Var}(S^2)=2m(m{-}1)+24q \equiv D,$$
$$\operatorname{Cov}(T,S)=6\tau\ \ (\textstyle\sum_{e\in E}A_e=3\tau),\qquad
\operatorname{Cov}(T,S^2)=\operatorname{Var}(T)=4P_1+16q\equiv W .$$
$\operatorname{Cov}(T,S^2)=\operatorname{Var}(T)$ is Lemma A in disguise: all of
$T$'s codegree form sits inside $S^2$, while $\operatorname{Cov}(R_4,S)=0$
(disjoint pairs cannot close a triangle). All verified as integer identities.

## Proposition (quadratic projection = the conditioning correction)

The best linear predictor of $T$ from $\{1,S\}$ has variance exactly
$\operatorname{Cov}(T,S)^2/\operatorname{Var}(S)=36\tau^2/m$ — the crude constant.
Projecting onto $\{1,S,S^2\}$ (Gram matrix from Lemma B) gives, for **any graph**,
$$\boxed{\ \operatorname{Var}(\mathbb E[T\mid c])\;\approx\;V_q
=\frac{36\tau^2 D-72\tau^2 W+m\,W^2}{m\,D-36\tau^2}\ },\qquad
W=4P_1+16q,\ D=2m(m{-}1)+24q,$$
an $O(n^3)$ closed form. The only approximation is that $\mathbb E[T\mid c]$ is
near-quadratic in $c$ — **measured**: $V_q$ captures 94–99% of the exact value on
every instance; a quartic captures 99%. For triangle-free graphs $V_q=W^2/D>0$ —
the previously unexplained "small residual" (random 3-regular:
$W^2/D\approx(12n)^2/4.5n^2\approx 32$).

**ER(n,p) leading order.** With $m\simeq pn^2/2$, $\tau\simeq p^3n^3/6$,
$q\simeq p^4n^4/8$, $W\simeq 2p^4n^4$, $D\simeq n^4(p^2/2+3p^4)$:
$$\operatorname{Var}(\mathbb E[T\mid c])\simeq \kappa(p)\,\frac{\tau^2}{m},\qquad
\kappa(p)=\frac{18\left(\tfrac12+2p-p^2\right)}{\tfrac14+\tfrac32 p^2-p^3},$$
with $\kappa(0.5)=45$, $\kappa(0.25)=51.4$, $\kappa\to36$ as $p\to0$ or $1$. The
measured "≈60" is the **finite-n value of $V_q$** (at n=12–14 the $4P_1$ term in
$W$ is not yet negligible), not a universal constant.

## Numerical verification — family means, 15 instances each (VB = exact Var(E[T|c]))

| family   | $V_q$/VB | deg-4 proj/VB | measured VB·m/τ² | closed-form $V_q$·m/τ² | crude |
|----------|---------|---------------|------------------|------------------------|-------|
| ER(0.5)  | 0.960   | 0.992         | 61.0             | 58.6                   | 36    |
| BA(k=4)  | 0.962   | 0.991         | 58.9             | 56.6                   | 36    |
| ER(0.25) | 0.984   | 0.996         | (τ tiny; ratio meaningless) | —           | —     |

Per instance $V_q$/VB ranges 0.94–0.997; the closed form reproduces the ≈60 to
~4% and explains its per-instance spread (measured ratios 45–84 track $V_q$).
For ER(0.25) the non-triangle $W^2/D$ part dominates, yet $V_q$ still captures
98% — one formula subsumes both regimes.

## The remainder (~4%) and its precise obstruction

The gap $VB-V_q$ is the degree-≥3 component of $\mathbb E[T\mid c]$. The cubic
moment is also closed-form (verified exactly on all instances):
$$\operatorname{Cov}(T,S^3)=6\tau(3m-2)+12\Phi,\qquad
\Phi=\sum_{j<k}A_{jk}\,P^{(3)}_{jk},$$
$P^{(3)}_{jk}$ = #simple 3-paths $j\to k$ ($\Phi$ ≈ pentagon-type counts). All
that resists closed form is the cubic/quartic **Gram matrix**: $\mathbb E[S^5]$,
$\mathbb E[S^6]$ need 5-/6-edge even-subgraph counts (5-, 6-cycles, thetas,
triangle pairs) — polynomial-time, so the degree-4 projection (99%) is in
principle closed-form too, but the extra 3% does not justify the bookkeeping.

## What this buys the paper

Since $V_2=\operatorname{Var}(T)-\operatorname{Var}(\mathbb E[T\mid c])$, at the
same fidelity $V_2\approx W-V_q$: the density law becomes a per-instance formula
from four subgraph counts $(m,\tau,q,P_1)$. For ER,
$V_2\simeq 2n^4p^4\big(1-p\,\kappa(p)/36\big)$, refining $2n^4p^4(1-p)$.

## Paper-ready LaTeX

```latex
\begin{lemma}[Quartic split]\label{lem:quarticsplit}
Let $S=\sum_{e\in E}X_e=m-2c$ with $X_{(j,k)}=s_js_k$. Then
$T=S^2+m-R_4$ with $R_4=\sum_{e\ne f,\,e\cap f=\emptyset}X_eX_f$ (ordered pairs).
Consequently $\Var_{S_v}[T]=\Var_{S_v}[R_4]$ for every cost class, and
$\Cov(T,S)=6\tau$, $\Cov(T,S^2)=\Var(T)=4P_1+16q$, $\Var(S)=m$,
$\E[S^3]=6\tau$, $\Var(S^2)=2m(m-1)+24q$, where $\tau,q,P_1$ count triangles,
$4$-cycles, and edge pairs sharing a vertex.
\end{lemma}

\begin{proposition}[Conditioning correction]\label{prop:condcorr}
Projecting $\E[T\mid c]$ onto quadratics in $c$ yields, with $W=4P_1+16q$ and
$D=2m(m-1)+24q$,
$\Var(\E[T\mid c])\approx(36\tau^2D-72\tau^2W+mW^2)/(mD-36\tau^2)$;
the linear part alone is $36\tau^2/m$. Numerically the quadratic projection
captures $94$--$99\%$ of the exact $\Var(\E[T\mid c])$ (E015), and for
$G(n,p)$ it equals $\kappa(p)\tau^2/m+o(n^4)$ with
$\kappa(p)=18(\tfrac12+2p-p^2)/(\tfrac14+\tfrac32p^2-p^3)$, e.g.\
$\kappa(0.5)=45$; the remainder is the degree-$\ge3$ component of
$\E[T\mid c]$, entering through $\Cov(T,S^3)=6\tau(3m-2)+12\sum_{j<k}A_{jk}P^{(3)}_{jk}$.
\end{proposition}
```
