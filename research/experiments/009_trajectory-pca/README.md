# 009 — Is the trajectory low-rank, and is the cost-class subspace aimed at it?

**Question:** For p = 20 ramp schedules, what fraction of the QAOA
trajectory's energy lies in the (m+1)-dimensional cost-class subspace
(E_cc), and how many trajectory-PCA dimensions (k_match) achieve the same —
i.e., when the proxy degrades, is the state genuinely high-dimensional, or
do few dimensions suffice but cost classes aren't them?

**Answer: "Wrong subspace," decisively: p=20 trajectories are extraordinarily
low-rank (k99 ≈ 2–4 for small/moderate ramps, ≤ 13 even at large), and a 2–3
dimensional trajectory-PCA subspace captures as much energy as the entire
(m+1 ≈ 22–46)-dimensional cost-class subspace — the proxy's subspace is
10–20× larger than the trajectory's effective dimension yet mis-aimed enough
to lose 4–21% of the energy outside the small-angle regime.**

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

Means over 10 instances at n=14 (n=12 similar; full table in `results.csv`):

| ramp     | E_cc range (across families) | k_match | k90     | k99      | m+1   |
|----------|------------------------------|---------|---------|----------|-------|
| small    | 0.974 – 0.989                | 2.0     | 1.0–1.7 | 2.1–3.2  | 22–46 |
| moderate | 0.941 – 0.966                | 2.1–2.5 | 2.0     | 3.0–4.1  | 22–46 |
| large    | 0.792 – 0.893                | 2.1–2.8 | 2.7–3.5 | 6.1–12.6 | 22–46 |

1. **QAOA ramp trajectories are tiny-rank objects**: ~99% of a 21-state,
   2¹⁴-dimensional trajectory fits in 2–4 PCA dimensions (slow ramps make
   consecutive states nearly parallel). Even at large angles, 13 dimensions
   suffice.
2. **The cost-class subspace is dimensionally lavish but mis-aimed**:
   k_match ≈ 2–3 everywhere — a 2–3-dim subspace does the job of all m+1
   cost-class dimensions. When the proxy degrades, it is *not* because the
   state needs more dimensions; the fixed cost-class frame simply drifts away
   from the trajectory's few directions.
3. Family ordering of E_cc again tracks density (ER(0.5)/BA(k=4) lowest),
   consistent with H-density and exp 003/007.
4. **Outlook this licenses (§7, careful not to overclaim):** instance-adapted
   low-rank compressions (e.g., subspaces built from a few cheap short-depth
   states) could in principle outperform the homogeneous proxy at equal
   dimension; the obstruction for parameter setting is that trajectory-PCA
   requires the very states the proxy is meant to avoid computing.

## Caveats

- Trajectory PCA is an *operationally unachievable* benchmark (it sees the
  answer before compressing); it is used to discriminate narratives, not as a
  competing method.
- E_cc measures state energy, not parameter-setting quality; exps 002/004
  cover the latter.
- Linear-ramp schedules only; non-ramp schedules with large per-layer angle
  changes would raise effective rank.

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/009_trajectory-pca/run.jl
```

Seed 20260611. Output: `results.csv`. Smoke test: prefix `E1_SMOKE=1`.
