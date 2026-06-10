# 002 — How much headroom does proxy parameter setting have above trivial baselines?

**Question:** Per graph family, how much approximation-ratio headroom exists
between trivial baselines (uniform random, random balanced partition) and the
best grid-achievable real QAOA — and how much of it does the exact-compression
proxy capture?

**Answer: PENDING — full run in progress (420 instances).**

## Why this matters

This is the metric-design experiment for the Phase-1 go/no-go gate. The
internship report observed that random balanced partitions already achieve
~75–85% AR on ER(0.5), so raw AR flatters any method. All downstream
experiments use **regret** = AR_ceiling − AR_proxy and **value-added** =
AR_proxy − AR_balanced instead. If value-added ≈ 0 within noise on *all*
families at p ≤ 3, parameter setting itself is ill-posed at this scale and the
paper reframes around state fidelity.

## Method

For each of 7 families (ER(0.5), ER(0.25), BA(k=2), BA(k=4), WS(k=4, β=0.1),
WS(k=4, β=0.5), 3-regular) × n ∈ {12, 14} × 30 instances:

- AR of uniform random measurement (mean cost / optimum) and the exact mean
  over all balanced partitions.
- True-QAOA grid ceiling: best AR over a 40×40 (γ, β) grid at p=1, and over an
  8⁴ linear-ramp endpoint grid (γ₁, γf, β₁, βf) at p=3.
- Proxy-set AR: real QAOA evaluated at the argmax of the **exact-compression
  proxy** (same-instance empirical N(c';d,c) and empirical P(c), per
  experiment 001 this is exactly the projected evolution — no model error) on
  the same grids.

The proxy here is the *best case* for the homogeneous heuristic: any analytical
or fitted N can only add model error on top.

## Result

*(pending; smoke run at n=10 suggests the proxy sits close to the p=1 ceiling
on all families, with a larger gap at p=3 ramps — do not cite, coarse grid)*

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/002_baselines-and-headroom/run.jl
```

Seed 20260611; per-instance seeds derived arithmetically (see script). Output:
`results.csv`. Smoke test: prefix `E1_SMOKE=1` (~1 min).
