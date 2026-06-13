# The V₂ density law: an exact anatomy (research note)

*Written 2026-06-13. Attacks the open problem flagged in the paper's Discussion:
"an analytic estimate of the variance functional V₂ per random-graph family,
which would turn the empirical density law into a theorem, is open and looks
tractable." Backed by experiment 015 (`research/experiments/015_v2-density-law/`),
which verifies the two exact lemmas to 1e-10 and the cubic law to <0.3% on 280
instances. This note gives the derivations and a paper-ready LaTeX block. Where
something is heuristic rather than proved, it says so.*

## Recap: what V₂ is

The cubic corollary (paper §4) writes the first-layer MaxCut leakage from
|+⟩^⊗n as
$$\lambda_1 = \frac{|\beta|\gamma^2}{8}\sqrt{V_2}\,(1+O(\gamma)+O(\beta)),\qquad
V_2 = 2^{-n}\sum_v M_v \operatorname{Var}_{S_v}[s_2],\quad s_2(y)=\sum_{i=1}^n c(y\oplus e_i)^2.$$
So V₂ is the within-cost-class variance of the sum of *squared neighbor costs*.
The paper measures √V₂ ∝ m ("density drives compression error") but treats V₂ as
a black box. We open it.

Notation: spins $s_j=(-1)^{y_j}\in\{\pm1\}$; an edge $(i,j)$ is cut iff
$s_is_j=-1$; $\deg(i)$, $\kappa_i(y)=$ #cut edges at $i$;
$\delta_i(y)=c(y\oplus e_i)-c(y)$; $A_{jk}=$ #common neighbors of $j,k$ (codegree),
with $A_{jj}=\deg(j)$; $\tau=$ #triangles.

## Lemma 1 (reduction to neighbor-spin sums)

**Claim.** $s_2(y) = (n-8)\,c(y)^2 + 4m\,c(y) + T(y)$ with
$T(y)=\sum_i \delta_i(y)^2$. Hence the first two terms are *class-constant* and
$$\boxed{\,V_2 = 2^{-n}\sum_v M_v \operatorname{Var}_{S_v}[T]\,.}$$

**Proof.** $c(y\oplus e_i)=c(y)+\delta_i$, so
$c(y\oplus e_i)^2=c(y)^2+2c(y)\delta_i+\delta_i^2$. Summing over $i$,
$s_2=n\,c^2+2c\sum_i\delta_i+T$. The neighbor-cost Lemma (paper Lemma 1)
gives $\sum_i c(y\oplus e_i)=(n-4)c+2m$, i.e. $\sum_i\delta_i=2m-4c$. Thus
$2c\sum_i\delta_i=4mc-8c^2$ and $s_2=(n-8)c^2+4mc+T$. The first two terms depend
on $y$ only through $c(y)$, so they are constant on each class $S_v$ and drop out
of the within-class variance. ∎

*Interpretation.* Only the **squared single-flip cost changes** $\delta_i^2$
matter — not the raw neighbor costs. This is the precise sense in which leakage
is a *local* (single-flip) second-moment effect.

## Lemma 2 (exact unconditional variance = sum of squared codegrees)

**Claim.** $T(y)-2m = 2\sum_{j<k}A_{jk}\,s_js_k$, and over **uniform** $y$,
$$\boxed{\;\mathbb E_y[T]=2m,\qquad \operatorname{Var}_y(T)=2\sum_{j\neq k}A_{jk}^2\;}$$
— an $O(n^2{+}m)$ closed form needing no $2^n$ enumeration.

**Proof.** $\kappa_i=\tfrac{\deg i}{2}-\tfrac12\sum_{j\sim i}s_is_j$, so
$\delta_i=\deg(i)-2\kappa_i=s_i\sum_{j\sim i}s_j$ and (using $s_i^2=1$)
$\delta_i^2=(\sum_{j\sim i}s_j)^2$. Then
$T=\sum_i\sum_{j,k\sim i}s_js_k=\sum_{j,k}A_{jk}s_js_k$ where
$A_{jk}=\#\{i:i\sim j,\,i\sim k\}$. The diagonal $j=k$ gives
$\sum_j A_{jj}s_j^2=\sum_j\deg(j)=2m$; the rest is
$2\sum_{j<k}A_{jk}s_js_k$. Over uniform spins $\mathbb E[s_js_k]=0$ ($j\neq k$) so
$\mathbb E[T]=2m$, and
$\operatorname{Var}(T)=4\sum_{j<k}\sum_{l<p}A_{jk}A_{lp}\,\mathbb
E[s_js_ks_ls_p]=4\sum_{j<k}A_{jk}^2$, since $\mathbb E[s_js_ks_ls_p]=1$ only when
$(j,k)=(l,p)$. ∎

Both lemmas are verified to relative tolerance $10^{-10}$ on all 280 instances
(exp 015 asserts them).

## The unconditional variance for Erdős–Rényi (exact leading order)

For $G(n,p)$, $A_{jk}\sim\mathrm{Binomial}(n-2,p^2)$ for $j\neq k$ (each other
vertex is a common neighbor independently w.p. $p^2$). Hence
$\mathbb E[A_{jk}^2]=(n-2)p^2(1-p^2)+((n-2)p^2)^2$ and
$$\mathbb E_{\text{graph}}\!\big[\operatorname{Var}_y(T)\big]
=2n(n-1)\big[(n-2)p^2(1-p^2)+((n-2)p^2)^2\big]
\;\xrightarrow{\;n\to\infty\;}\;2n^4p^4 = 8p^2m^2,$$
using $m\approx pn^2/2$. So the *unconditional* $\sqrt{\operatorname{Var}(T)}/m
\to \sqrt8\,p\approx 2.83p$ — it grows **linearly in the density** $p$, i.e.
super-linearly in $m$ at fixed $n$. (Measured: ER(0.5) ≈ 1.4, ER(0.25) ≈ 0.87 at
$n=16$ — consistent with finite-$n$ corrections.)

## Why √V₂ ∝ m: the conditioning cancellation

V₂ is the *within-class* variance, a fraction
$\rho_{\rm cond}=V_2/\operatorname{Var}_y(T)$ of the unconditional one. Writing
$Q=\sum_{(ij)\in E}s_is_j=m-2c$, conditioning on the cut value $c$ fixes $Q$.
Only pairs $(j,k)$ that are **edges and share a neighbor** (i.e. triangles) have
their $\mathbb E[s_js_k\mid c]$ pulled away from 0; non-edge codegree pairs stay
$\approx0$. Since $\sum_{(j,k)\in E}A_{jk}=3\tau$ (each triangle contributes one
common neighbor to each of its three edges),
$$\mathbb E[T\mid c]\approx 2m + \frac{6\tau}{m}(m-2c),\qquad
\operatorname{Var}(\mathbb E[T\mid c])\approx\Big(\frac{12\tau}{m}\Big)^2\!\operatorname{Var}(c)
\approx \frac{36\tau^2}{m},$$
using $\operatorname{Var}_y(c)\approx m/4$ (edge spin-products are pairwise
uncorrelated). Therefore
$$V_2 \approx \operatorname{Var}_y(T)-\frac{36\tau^2}{m}
= 2\sum_{j\neq k}A_{jk}^2 - \frac{36\tau^2}{m}.$$
For ER, $\tau\approx\binom{n}{3}p^3\approx n^3p^3/6$, giving
$36\tau^2/m\approx 2n^4p^5$, so
$$V_2 \approx 2n^4p^4(1-p),\qquad \frac{\sqrt{V_2}}{m}\approx 2\sqrt2\,p\sqrt{1-p}.$$
**The factor $\sqrt{1-p}$ is the point.** The unconditional $\propto p$ growth is
partly cancelled by the triangle-driven suppression $\sqrt{1-p}$, which is
*stronger for denser graphs*. The result is far flatter in density than the
unconditional variance — the empirical density law.

**Honest status of the constant.** The exact lemmas (1, 2) and the ER
$\mathbb E[\operatorname{Var}_y(T)]=8p^2m^2$ are rigorous. The conditioning
estimate $\rho_{\rm cond}\approx 1-p$ captures the **direction** (denser ⇒ more
suppression) but **overshoots**: measured $\rho_{\rm cond}$ at $n=16$ is ≈0.37
(ER(0.5)) and ≈0.60 (ER(0.25)) vs the heuristic 0.5 / 0.75. The crude step is
$\mathbb E[s_js_k\mid c]\approx Q/m$ (exchangeable-edge approximation) plus
dropping induced non-edge correlations. Pinning this constant — equivalently,
the exact $\operatorname{Var}(\mathbb E[T\mid c])$ — is the remaining open piece.
Experiment 015 confirms the **cancellation** numerically: across the seven
families the unconditional $\sqrt{\operatorname{Var}(T)}/m$ spans ~2.0–2.35×,
while $\sqrt{V_2}/m$ spans only ~1.35–1.58× (n=12–18).

## What this buys the paper

1. **The black box is opened.** Leakage's driver is the within-class variance of
   $T=\sum_i(\sum_{j\sim i}s_j)^2$, whose unconditional value is *exactly*
   $2\sum_{j\neq k}A_{jk}^2$ — squared codegrees. "Density drives compression
   error" sharpens to "squared codegrees, conditioned on cost, drive it."
2. **The density law gets a mechanism**, not just a fit: super-linear codegree
   growth × triangle-driven conditioning suppression ≈ linear in $m$.
3. **The open problem shrinks** to one quantity, $\operatorname{Var}(\mathbb
   E[T\mid c])$ (the triangle-mediated conditioning correction).

## Paper-ready LaTeX (proposed addition after Corollary 1 / cubic law)

> Insert as two lemmas in §4 (theory) or fold into the Discussion's open-problem
> paragraph. They are as rigorous and machine-verified as the existing Lemma.

```latex
\begin{lemma}[Variance reduction]\label{lem:v2reduce}
For MaxCut, $s_2(y)=(n-8)c(y)^2+4m\,c(y)+T(y)$ with
$T(y)=\sum_{i}(c(y\oplus e_i)-c(y))^2$, so the first two terms are class-constant
and $V_2 = 2^{-n}\sum_v M_v \Var_{S_v}[T]$.
\end{lemma}

\begin{lemma}[Codegree form]\label{lem:codegree}
Writing $s_j=(-1)^{y_j}$ and $A_{jk}$ for the number of common neighbors of $j,k$,
$T(y)-2m=2\sum_{j<k}A_{jk}\,s_js_k$; hence over the uniform distribution
$\mathbb E[T]=2m$ and $\Var(T)=2\sum_{j\neq k}A_{jk}^2$. For $G(n,p)$,
$\mathbb E[\Var(T)]=2n(n-1)[(n-2)p^2(1-p^2)+((n-2)p^2)^2]\to 8p^2m^2$.
\end{lemma}

\begin{remark}
$V_2$ is the within-class (cost-conditioned) version of $\Var(T)$. The
unconditional $\sqrt{\Var(T)}/m\to\sqrt8\,p$ grows with density, but conditioning
on $c$ removes a triangle-driven fraction that grows with density too; the two
nearly cancel, leaving $\sqrt{V_2}\propto m$ across families (E015). An exact
$\Var(\mathbb E[T\mid c])$ would close the density law.
\end{remark}
```
