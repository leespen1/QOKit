# E044 — Does the composed recipe (sampled N + normalized) survive n=16?

**Question**: Does E041's headline (sampled S=10 + normalized objective is
the best proxy combination at both depths, beating even exact+normalized)
replicate one size up, at n=16, including the S=3 cell E042 found strongest?

**Answer**: **PENDING** (full run in progress).

## Method

7 families × n=16 × 5 instances (E041 seed scheme). Variants
{exact, samp10, samp3} × {raw, normalized}; p=1 grid (40×40) and p=3
linear-ramp grid (8⁴); regret vs the true statevector grid ceiling.
Mirrors E040's stability role for the earlier findings.

## Result

PENDING.

## Caveats

- 5 instances per family (35 total), same as E040's n=16 budget.
- Single sampling draw per (instance, S); E042 quantified draw variance.

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/044_composed-recipe-n16/run.jl   # ~40 min
E44_SMOKE=1 julia --project research/experiments/044_composed-recipe-n16/run.jl              # quick check
```
