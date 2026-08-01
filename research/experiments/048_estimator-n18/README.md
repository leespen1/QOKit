# E048 — The n=18 point: does the estimator's depth deficit keep growing?

**Question**: E047's depth-split verdict (p=1 advantage grows with n, p=3
deficit grows with n) rests on two sizes; does n=18 confirm the trends?

**Answer**: **PENDING** (full run in progress, ~2-3 h).

## Method

7 families × n=18 × 3 instances (E041 seed scheme), exact / samp10 /
poly(S=10, K=200), both objectives, p=1 grid (40×40) and p=3 linear-ramp
grid (8⁴), regret vs the statevector grid ceiling (2^18). Harness derived
from E047 part B; estimator.jl include()d from E046 unmodified.

## Result

PENDING.

## Caveats

- 3 instances per family (21 total): trend-check budget, not a
  full-precision cell.

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/048_estimator-n18/run.jl   # ~2-3 h
E48_SMOKE=1 julia --project research/experiments/048_estimator-n18/run.jl              # quick check
```
