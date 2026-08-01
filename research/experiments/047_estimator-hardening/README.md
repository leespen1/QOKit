# E047 — Hardening the polynomial-time N estimator

**Question**: Is E046's estimator robust to its hyperparameters (S members,
K profile samples), does its p=1 advantage persist at n=16, and how does
its wall-clock scale to n=24?

**Answer**: **Split by depth. At p=1 the estimator is robust and improves
with size (n=16: mean regret 0.0084, better than the enumeration sampler
on 29/35; at n=12, LESS profile sampling is better, K=50 giving 0.0061 —
the noise-regularization mechanism again). At p=3 the gap to the
enumeration sampler does not close under any (S,K) tested and WIDENS with
n (~0.005 pooled at n≤14 → +0.022 at n=16, poly better only 7/35). Timing
scales excellently: n=24 estimation in ~1.3 s, ~2.6M Wang-Landau steps,
zero missed classes everywhere.**

## Method

Part A: n=12, 7 families × 5 instances, poly at (S,K) ∈
{(3,50),(10,50),(10,200),(10,800),(30,200)} vs exact and samp10, both
depths, both objectives ((10,200) reproduces E046's rng draws exactly).
Part B: same harness at n=16 with poly(10,200). Part C: timing at
n ∈ {16,20,24} × 3 ER(0.5) instances with instrumented Wang-Landau
step/stage counts; enumeration timed for n≤20 only. Fixed seeds
(SEED=20260611); estimator.jl include()d from E046 unmodified.

## Result

Part A (n=12, mean regret_norm): p=1 — K=50 best (0.0061 at S=10, 0.0120
at S=3), K=200 0.0142, K=800 0.0242, S=30 0.0257, vs samp10 0.0187 and
exact 0.0435: more sampling (less noise) hurts p=1, exactly as E042
predicts. p=3 — all poly variants cluster 0.027–0.030 vs samp10 0.0224;
K=800 closes almost nothing (0.0274); S=3+K=50 degrades (0.0437).

Part B (n=16): p=1 poly 0.0084 (max 0.0198) vs samp10 0.0183, exact
0.0609 — the advantage grows with n. p=3 poly 0.0453 vs samp10 0.0238 —
the deficit also grows with n (paired +0.0216, poly better 7/35).
Per-family pooled: poly wins ER(0.5), BA(k=2), BA(k=4); loses the rest,
driven by p=3.

Part C: t_poly ≈ 1.2–1.5 s at n=24 (m ≈ 119–147), ~2.4–2.6M WL steps, 10
stages; enumeration already slower at n=20 (E046) and infeasible-in-kind
beyond. classes_missed = 0 on every instance at every size.

## Caveats

- 5 instances per family; single draw per (instance, S, K).
- n=16 is one size beyond validation; the depth-deficit trend
  (0.005 → 0.022) rests on two points and merits an n=18 check before any
  stronger claim.
- No statevector truth at n=24 (timing only), by design.

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/047_estimator-hardening/run.jl   # ~50 min
E47_SMOKE=1 julia --project research/experiments/047_estimator-hardening/run.jl              # quick check
```
