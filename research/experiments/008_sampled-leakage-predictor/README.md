# 008 — Can leakage be estimated from a few sampled bitstrings per cost class?

**Question:** How accurately does Theorem 3's variance formula, evaluated on
S ∈ {2, 5, 10, 25} stratified samples per cost class, reproduce the exact
per-class leakage λ(c′) and its aggregate η_F — i.e., can the paper's central
diagnostic be computed at scales where the O(4ⁿ) exact distribution is
unreachable?

**Answer: PENDING — full run in progress. (Smoke run already established one
planned result: at S = full enumeration the variance formula reproduces the
statevector-computed leakage to 10⁻¹⁶, verifying Theorem 3's identity via two
independent computational routes.)**

## Why this matters

- Re-verifies old hypothesis H-sample ("S = 5–10 suffices") in the quantity
  that matters for this paper — leakage — rather than entrywise N estimation.
- A sampled leakage estimator is the enabler for Phase-3 scale-up (n = 18–20),
  where exact homodist/leakage computation is expensive.
- The S = 0 (full enumeration) rows at n ≤ 12 are a machine-precision check of
  Theorem 3 itself (λ(c′)² as class-size-weighted within-class variance of the
  f_d(β)-weighted neighborhood profile g).

## Method

7 families × n ∈ {12, 14} × 10 instances (shared seeds), 3 (γ, β) points
(small → moderate). Exact λ(c′) per attained class via one statevector layer +
projection; sampled λ̂(c′) via per-class stratified sampling of g_{c′}(y) (one
O(2ⁿ) neighborhood-histogram scan per sampled y), with the proper population/
sample variance conventions. Reported: relative error of η_F and median
relative error of per-class λ.

## Result

*(pending)*

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/008_sampled-leakage-predictor/run.jl
```

Seed 20260611. Output: `results.csv`. Smoke test: prefix `E1_SMOKE=1`.
