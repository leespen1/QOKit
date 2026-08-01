# E048 — The n=18 point: does the estimator's depth deficit keep growing?

**Question**: E047's depth-split verdict (p=1 advantage grows with n, p=3
deficit grows with n) rests on two sizes; does n=18 confirm the trends?

**Answer**: **The depth deficit is real and persists (p=3: poly trails
samp10 by +0.026 paired, better only 2/21 — three sizes now: 0.005 →
0.022 → 0.026), but the p=1 advantage does NOT keep growing: at n=18 poly
and samp10 tie (0.0154 vs 0.0143, 11w/10l). Honest summary: the estimator
matches the enumeration sampler at p=1 across n=12-18 and trails it at
depth. The composed recipe itself keeps beating exact at n=18 (samp10+norm
better than exact+norm on 20/21 at p=1, 17/19 decided at p=3). Zero missed
classes.**

## Method

7 families × n=18 × 3 instances (E041 seed scheme), exact / samp10 /
poly(S=10, K=200), both objectives, p=1 grid (40×40) and p=3 linear-ramp
grid (8⁴), regret vs the statevector grid ceiling (2^18). Harness derived
from E047 part B; estimator.jl include()d from E046 unmodified.

## Result

Mean regret_norm over 21 instances (raw in parentheses):

| variant | p=1 | p=3 |
|---|---|---|
| exact | 0.0662 (0.043) | 0.0691 (0.089) |
| samp10 | 0.0143 (0.041) | 0.0268 (0.082) |
| poly S=10 K=200 | 0.0154 (0.184) | 0.0528 (0.296) |

Paired poly−samp10: p=1 +0.0011 (11w/10l), p=3 +0.0260 (2w/19l).
Paired samp10−exact: p=1 −0.0519 (20w/1l), p=3 −0.0423 (17w/2l) — the
E041/E044 composition result extends to n=18. The p=1 size trend for poly
(n=12: −0.004, n=16: −0.010, n=18: +0.001 vs samp10) is not monotone;
"matches" is the defensible claim. classes_missed = 0 everywhere.

## Caveats

- 3 instances per family (21 total): trend-check budget, not a
  full-precision cell.

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/048_estimator-n18/run.jl   # ~2-3 h
E48_SMOKE=1 julia --project research/experiments/048_estimator-n18/run.jl              # quick check
```
