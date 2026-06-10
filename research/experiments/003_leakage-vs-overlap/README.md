# 003 — Does accumulated leakage track the true-vs-proxy overlap deficit?

**Question:** For deep (p = 20) linear-ramp schedules, how tightly does the
accumulated per-layer leakage Σλ_ℓ bound the actual distance ‖ψ_p − φ_p‖
between the true and compressed trajectories (Theorem 2), across graph
families and ramp magnitudes?

**Answer: PENDING — full run in progress (840 trajectories).**

## Why this matters

Theorem 2 (‖ψ_p − φ_p‖ ≤ Σ_ℓ λ_ℓ) is the paper's quantitative link between
per-layer leakage and end-to-end proxy fidelity. If the bound is loose by >10×
systematically, it is true but vacuous and the paper leans on measured λ_ℓ
directly; if it is within a few ×, it justifies using cheap per-layer leakage
as the central diagnostic.

## Method

7 families × n ∈ {12, 14, 16} × 10 instances × 4 linear-ramp schedules
(small/moderate/large/extreme angle magnitudes), p = 20. For every layer we
record λ_ℓ, Σλ, ‖ψ_ℓ − φ_ℓ‖, |⟨ψ_ℓ|φ_ℓ⟩|, and ‖φ_ℓ‖ via
`compressed_qaoa_trajectory` (instance seeds shared with experiment 002).

## Result

*(pending; smoke run at n=10: bound holds everywhere, slack ~3–5× at p=20,
overlap ranges 0.95+ for small ramps down to ~0 for extreme ones — exactly the
dynamic range we want. Do not cite.)*

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/003_leakage-vs-overlap/run.jl
```

Seed 20260612. Output: `results.csv` (long format, one row per layer).
Smoke test: prefix `E1_SMOKE=1`.
