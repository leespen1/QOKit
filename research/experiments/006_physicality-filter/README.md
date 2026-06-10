# 006 — Does a physicality cap (predicted ⟨C⟩ ≤ m) repair the analytical proxy without collateral damage?

**Question:** Restricting the analytical PaperProxy's argmax to grid points
whose predicted ⟨C⟩ does not exceed the number of edges (an impossibility
bound, unlike experiment 005's calibration-sensitive norm tolerance) — does
that fix the dense-ER(0.5) spurious-peak failure while leaving the
already-good sparse-family argmaxes untouched?

**Answer: PENDING — full run queued.**

## Why this matters

Experiment 005 established that the analytical model's *calibration* (values,
norms) is broken on sparse families even where its *argmax location* is nearly
perfect, so calibration-based filters are catastrophic. A cap at the physical
maximum rejects only predictions that are certainly wrong, never a well-located
peak with merely inflated height — unless the inflation pushes it past m, which
is exactly the ER(0.5) failure mode. If this works, the paper's practical
recommendation becomes: *use the analytical N anywhere, but never trust a
prediction that exceeds the physical maximum.* A tighter heuristic cap
(0.75·m) is included as a sensitivity check; when a cap rejects the whole grid
the selector falls back to the raw argmax (recorded via frac columns).

## Method

Same 420 instances and grids as experiments 004/005. PaperProxy arm only:
AR at argmax under no cap / cap ⟨C⟩ ≤ m / cap ⟨C⟩ ≤ 0.75·m, at p=1 and p=3.

## Result

*(pending)*

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/006_physicality-filter/run.jl
```

Seed 20260611 (same instance set as 002/004/005). Output: `results.csv`.
Smoke test: prefix `E1_SMOKE=1`.
