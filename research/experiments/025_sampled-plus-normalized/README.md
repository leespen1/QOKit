# E025 — Does the recommended recipe compose? Sampled N + normalized objective

**Question**: Does the paper's recommended recipe — sampled N (S=10 per cost
class) *combined with* the normalized objective at depth — actually deliver
the regret its two separately-validated ingredients (E010: sampled N,
unnormalized; E020/E024: normalization, exact N) promise?

**Answer**: **PENDING** (full run in progress; a 123/140-instance partial run
from the 2026-07-22 loop was killed at the window boundary and is being rerun).

## Method

Standard instance set: 7 families × n ∈ {12, 14} × 10 instances, fixed seeds
(`instance_seed(fi, n, inst)`, base seed 20260611). For each instance and each
depth (p=1 grid; p=3 linear-ramp grid), compute the true QAOA ceiling by
statevector simulation, then the regret of the grid point selected by each of
the four proxy objectives: {exact N, sampled N (S=10)} × {unnormalized,
normalized}. Sampling uses an independent RNG stream (seed + 777).

## Result

PENDING.

## Caveats

- Grid ceilings, not continuous optima (E018 showed these are tight).
- S=10 per class is the single sampling budget tested; no S-sweep here.

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/025_sampled-plus-normalized/run.jl   # ~15 min
E25_SMOKE=1 julia --project research/experiments/025_sampled-plus-normalized/run.jl              # quick check
```
