# 015 — Why does compression error scale with density? An exact anatomy of V₂

**Question:** The cubic corollary says first-layer MaxCut leakage is
λ₁ = (|β|γ²/8)·√V₂, and the paper measures √V₂ ∝ m ("density, not ER-ness")
but leaves an analytic estimate of V₂ open. Can we pin down what V₂ *is*,
and explain mechanistically why it tracks edge count?

**Answer: Yes. Two exact reductions (verified to 1e-10 on 280 instances) cut V₂
down to the within-cost-class variance of T(y)=Σ_i(Σ_{j∼i}s_j)², whose
*unconditional* variance is exactly 2·Σ_{j≠k}A_{jk}² (A = codegree, #common
neighbors). That unconditional variance grows super-linearly with density
(∝ density²·m² for Erdős–Rényi), but conditioning on cost removes a
density-dependent fraction ρ_cond (≈0.35 for dense ER(0.5), ≈0.80 for sparse
3-regular), and the two effects cancel: the family spread of √leakage/m drops
from ~2× (unconditional) to ~1.4× (the actual V₂) — the density law is this
near-cancellation, not a coincidence. The cubic law itself is essentially exact:
λ₁/λ₁_pred ∈ [0.997, 1.002] over every instance.**

## Why this matters

The paper's §4 mechanism (compression error is third-order, ∝ βγ²·√V₂) and its
§5 headline ("density drives compression error") both rest on V₂, which was a
black box. This experiment opens it: it gives the *exact* graph functional
behind leakage and a quantitative reason the density law holds across families.
It directly advances the open problem flagged in the Discussion ("an analytic
estimate of V₂ per family would turn the empirical density law into a theorem").

## The exact results (proved + verified)

Throughout, s = (−1)^{y} are ±1 spins and δ_i(y) = c(y⊕e_i) − c(y) =
deg(i) − 2κ_i(y) = Σ_{j∼i} s_j (κ_i = cut edges at i).

- **Lemma L1 (reduction).** s₂(y) = Σ_i c(y⊕e_i)² = (n−8)c(y)² + 4m·c(y) + T(y),
  where T(y) = Σ_i δ_i². The first two terms are class-constant (they depend on y
  only through c(y)), so they vanish under the within-class variance:
  **V₂ = 2⁻ⁿ Σ_v M_v Var_{S_v}[T].** Only the squared neighbor-spin sums matter,
  not the raw neighbor costs. (Uses the committed Lemma Σ_i c(y⊕e_i)=(n−4)c+2m.)
- **Lemma L2 (codegree closed form).** T(y) − 2m = 2·Σ_{j<k} A_{jk} s_j s_k with
  A_{jk} = #common neighbors of j,k, so over uniform y,
  **Var_y(T) = 2·Σ_{j≠k} A_{jk}².** (E[T] = 2m.) This is an O(n²)+O(m) closed form
  — no 2ⁿ enumeration — for the *unconditional* variance.
- **Cubic law C.** λ₁ = (|β|γ²/8)√V₂ to leading order. Measured at γ=β=0.05:
  λ₁/λ₁_pred ∈ [0.997, 1.002] across all 7 families, n=12–18 (mean 0.9985).

Both L1 and L2 are verified to relative tolerance 1e-10 on every one of the 280
instances (the run asserts them; it fails fast otherwise).

## The mechanism (measured)

V₂ is the *within-class* (conditional) variance, a fraction
ρ_cond = V₂ / Var_y(T) of the unconditional codegree variance:

| | dense ←——————————————→ sparse | |
|---|---|---|
| family | ER(0.5) · BA(4) · WS(0.5) · BA(2) · WS(0.1) · ER(0.25) · 3-reg | |
| ρ_cond (n=16) | 0.37 · 0.41 · 0.63 · 0.56 · 0.53 · 0.60 · 0.78 | denser ⇒ less survives |

The unconditional √Var(T)/m varies ~2.0–2.35× across families (it grows like
density·m for ER, super-linear in m); ρ_cond shrinks with density (more
triangles ⇒ the cut value c explains more of T); the product √(ρ_cond·Var(T))/m
= √V₂/m varies only ~1.35–1.58×. So conditioning roughly **halves** the family
spread, flattening leakage toward √V₂ ≈ const·m.

For Erdős–Rényi G(n,p) the unconditional piece is analytic: A_{jk} ~
Binomial(n−2, p²), so E[Var_y(T)] = 2·n(n−1)·[(n−2)p²(1−p²) + ((n−2)p²)²]
≈ 2n⁴p⁴ = 8p²·m² (leading order). The remaining piece — E[Var(E[T|c])], the
triangle-mediated conditioning correction — is the only ingredient still needed
for a closed-form E[V₂]; this experiment isolates and quantifies it.

## Results

`results.csv`: 280 rows (7 families × n∈{12,14,16,18} × 10 instances, the exp
002/004 instance set). Columns: family, n, inst, seed, m, V2, V2_from_s2,
varT_uncond, varT_codegree, SA2, rho_cond, lambda1, lambda1_pred, lambda_ratio.
`v2_anatomy.png` (from `make_figure.jl`): (a) measured λ₁ vs (βγ²/8)√V₂ — the
cubic law on the diagonal; (b) √Var(T)/m vs √V₂/m by family — the flattening.

## Caveats

- L1, L2, and the βγ²√V₂ law are **exact / leading-order theorems** (any graph);
  the density law √V₂ ∝ m and the ρ_cond trend are **measured** over the tested
  families and n=12–18. A closed-form E[V₂] is not yet derived (the conditioning
  correction is open, now isolated).
- V₂ here is the *first-layer* leakage coefficient from |+⟩ⁿ. Deeper-layer
  leakage uses the same variance identity (Thm 3) about a general class state;
  this experiment does not re-derive the constant there.
- Brute-force V₂ is O(n·2ⁿ); n≤18 on CPU. The codegree form (L2) is cheap but
  gives only the unconditional Var(T).

## Reproduce

```
julia --project research/experiments/015_v2-density-law/run.jl
julia --project research/experiments/015_v2-density-law/make_figure.jl
```

Instances seeded by `20260611 + 10000·fam_idx + 100·n + inst` (identical to
exp 002/004). Smoke test: `E15_SMOKE=1 julia --project .../run.jl`.
