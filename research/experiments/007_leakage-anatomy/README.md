# 007 — Leakage anatomy over the (γ, β) plane

**Question:** How does one-layer leakage out of the cost-class subspace vary
with (γ, β) per graph family — and in particular, is the proxy's trust region
predictable from leakage alone (does the analytical proxy's dense-ER spurious
peak sit in a high-leakage region)?

**Answer: PENDING — full run in progress.**

## Why this matters

Three threads converge here:
1. The filter arc (exps 004–006) showed value/norm-based guards cannot repair
   the analytical proxy; if high *leakage* marks the regions where any
   homogeneous proxy is untrustworthy, leakage itself defines the guard —
   and it is model-independent and computable from instance data.
2. H-density (exp 003): leakage tracked edge count at fixed angles; a (γ, β)
   map lets us test whether rescaling γ by m absorbs the density effect.
3. Theorem 3 predicts the small-angle scaling of leakage; this sweep provides
   the data to check the perturbative regime, and is the paper's planned
   "anatomy of leakage" figure (Fig. 1 in the skeleton).

## Method

7 families × n ∈ {12, 14} × 10 instances (seeds shared with exps 002/004–006),
24×24 grid over γ ∈ (0, π], β ∈ (0, π/2]. Per grid point: λ_uniform (one-layer
leakage from |+⟩ⁿ) and η_F = sqrt(Σ_{c'} λ(c')²) over all attained cost-class
states (model-independent operator-leakage profile). All exact statevector
computation via `apply_phase_gate!`/`apply_x_mixer!` + projection; no
distribution arrays.

## Result

*(pending)*

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/007_leakage-anatomy/run.jl
```

Seed 20260611. Output: `results.csv` (long format, one row per grid point).
Smoke test: prefix `E1_SMOKE=1`.
