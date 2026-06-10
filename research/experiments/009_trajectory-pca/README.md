# 009 — Is the trajectory low-rank, and is the cost-class subspace aimed at it?

**Question:** For p = 20 ramp schedules, what fraction of the QAOA
trajectory's energy lies in the (m+1)-dimensional cost-class subspace
(E_cc), and how many trajectory-PCA dimensions (k_match) achieve the same —
i.e., when the proxy degrades, is the state genuinely high-dimensional, or
do few dimensions suffice but cost classes aren't them?

**Answer: PENDING — full run in progress.**

## Why this matters

This experiment discriminates between two §4 narratives for proxy failure at
depth/large angles: "not low-rank" (no small subspace can hold the state — the
compression idea itself is doomed there) vs "wrong subspace" (a handful of
directions would suffice, but the cost-class directions are mis-aimed — the
compression idea survives with a better subspace, e.g. instance-adapted).
Pre-registered reading: if k_match ≪ m+1 while E_cc < 1, it's "wrong subspace."

## Method

7 families × n ∈ {12, 14} × 10 instances (shared seeds) × 3 ramps
(small/moderate/large), p = 20. Intermediate states via
`qaoa_statevector(...; return_intermediates=true)`; E_cc from
`project_onto_cost_classes`; optimal rank-k captured energy from the thin-SVD
spectrum of the stacked 2ⁿ × 21 trajectory matrix; k90/k99 = PCA dimensions
for 90%/99% of trajectory energy.

## Result

*(pending; smoke run at n=10 suggests trajectories are extremely low-rank —
k99 ≈ 3–5 — with E_cc ≈ 0.95–0.98 and k_match ≈ 2–3. Do not cite.)*

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/009_trajectory-pca/run.jl
```

Seed 20260611. Output: `results.csv`. Smoke test: prefix `E1_SMOKE=1`.
