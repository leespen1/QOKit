# 002 — How much headroom does proxy parameter setting have above trivial baselines?

**Question:** Per graph family, how much approximation-ratio headroom exists
between trivial baselines (uniform random, random balanced partition) and the
best grid-achievable real QAOA — and how much of it does the exact-compression
proxy capture?

**Answer: Real headroom exists on every family and the exact-compression proxy
captures most of it at p=1 (mean regret ≈ 0.03, value-added positive for every
single instance), while regret grows ~2–3× at p=3 ramps — so parameter setting
is well-posed, and depth is where the proxy starts leaving value behind.**

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

Means over 30 instances per cell (`results.csv` has per-instance rows).
VA = value-added = AR_proxy − AR_balanced; regret = AR_ceiling − AR_proxy.

| family        | n  | AR_balanced | VA (p=1) | min VA (p=1) | regret (p=1) | VA (p=3) | regret (p=3) |
|---------------|----|-------------|----------|--------------|--------------|----------|--------------|
| ER(0.5)       | 12 | 0.750       | 0.039    | 0.010        | 0.028        | 0.077    | 0.080        |
| ER(0.5)       | 14 | 0.764       | 0.031    | 0.015        | 0.030        | 0.064    | 0.080        |
| ER(0.25)      | 12 | 0.642       | 0.112    | 0.067        | 0.031        | 0.168    | 0.082        |
| ER(0.25)      | 14 | 0.646       | 0.098    | 0.068        | 0.036        | 0.150    | 0.089        |
| BA(k=2)       | 12 | 0.670       | 0.088    | 0.067        | 0.030        | 0.135    | 0.085        |
| BA(k=2)       | 14 | 0.662       | 0.086    | 0.077        | 0.038        | 0.142    | 0.080        |
| BA(k=4)       | 12 | 0.749       | 0.032    | 0.028        | 0.028        | 0.078    | 0.069        |
| BA(k=4)       | 14 | 0.746       | 0.042    | 0.037        | 0.028        | 0.070    | 0.087        |
| WS(k=4;b=0.1) | 12 | 0.719       | 0.076    | 0.071        | 0.028        | 0.123    | 0.063        |
| WS(k=4;b=0.1) | 14 | 0.703       | 0.082    | 0.075        | 0.031        | 0.112    | 0.084        |
| WS(k=4;b=0.5) | 12 | 0.693       | 0.076    | 0.068        | 0.031        | 0.129    | 0.077        |
| WS(k=4;b=0.5) | 14 | 0.679       | 0.084    | 0.068        | 0.035        | 0.120    | 0.092        |
| 3-regular     | 12 | 0.629       | 0.115    | 0.103        | 0.037        | 0.173    | 0.082        |
| 3-regular     | 14 | 0.629       | 0.121    | 0.104        | 0.038        | 0.165    | 0.092        |

Readings:

1. **Parameter setting is well-posed at this scale**: value-added is positive
   for all 420 instances individually (per-instance minimum 0.010, on dense
   ER(0.5) where the balanced floor is highest). The E1.1 falsification
   condition ("VA ≈ 0 within noise everywhere") is decisively not met.
2. **The balanced-partition floor confirms hypothesis H-baseline**: ~0.75–0.76
   on ER(0.5) and BA(k=4) (dense), leaving those families the least headroom;
   sparse families (ER(0.25), BA(k=2), 3-regular) have 2–3× more.
3. **At p=1 the exact-compression proxy is near-ceiling on *every* family**,
   including BA and WS (regret ≈ 0.03 across the board). The old log's
   "proxies fail on non-ER graphs" (H-ER-specific) therefore cannot be a
   compression-error story at p=1 — it must come from *model error* in the
   analytical/fitted N. To be tested directly in E1.3 by running PaperProxy N
   on the same instances.
4. **Regret grows to ~0.06–0.09 at p=3 ramps** on all families — consistent
   with leakage accumulating with depth (Theorem 2); E1.2 measures this
   directly.

## Caveats

- The p=3 "ceiling" is the best of the same 8⁴ ramp-endpoint grid the proxy
  chooses from (fair comparison, but a finer grid or free 2p-parameter
  optimization would raise both numbers).
- n ≤ 14 only; scale dependence is Phase 3's E3.1.
- The proxy here is the *exact-compression* (instance-empirical N and P) best
  case; analytical-N regret will be at least this large.

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/002_baselines-and-headroom/run.jl
```

Seed 20260611; per-instance seeds derived arithmetically (see script). Output:
`results.csv`. Smoke test: prefix `E1_SMOKE=1` (~1 min).
