# 008 — Can leakage be estimated from a few sampled bitstrings per cost class?

**Question:** How accurately does Theorem 3's variance formula, evaluated on
S ∈ {2, 5, 10, 25} stratified samples per cost class, reproduce the exact
per-class leakage λ(c′) and its aggregate η_F — i.e., can the paper's central
diagnostic be computed at scales where the O(4ⁿ) exact distribution is
unreachable?

**Answer: Yes — S = 5 samples per class estimate the aggregate leakage η_F to
a median 3.2% (p90: 8.9%) relative error, uniformly across families and
angles; family-level leakage differences are tens of percent, so S = 5–10
amply suffices for the paper's ranking claims (H-sample re-verified in the
quantity that matters). The S = full rows verify Theorem 3's identity to
4.5×10⁻¹⁵ over 210 instance/angle combinations.**

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

Relative error of η_F (pooled over 7 families × n ∈ {12,14} × 10 instances ×
3 angle points):

| S   | median   | p90    | max    |
|-----|----------|--------|--------|
| 2   | 5.6%     | 14.8%  | 33%    |
| 5   | 3.2%     | 8.9%   | 22%    |
| 10  | 2.2%     | 5.8%   | 13%    |
| 25  | 1.3%     | 3.7%   | 8%     |
| all | 5×10⁻¹⁶  | —      | 4×10⁻¹⁵ |

Per-class λ(c′) median relative error: 15% (S=2) → 7.8% (S=5) → 3.2% (S=25).
Accuracy is uniform across families (median at S=5: 2.2–5.0%) and across the
three (γ, β) points.

Consequences:
1. **Phase-3 enabler confirmed**: leakage at n = 18–22 can be estimated by
   S·|C| neighborhood scans (O(S·m·2ⁿ)) instead of any O(4ⁿ) computation,
   with errors far below the between-family differences being ranked.
2. **Theorem 3 is established numerically** (max 4.5×10⁻¹⁵ over 210 full-
   enumeration rows at n=12): the statevector route (evolve + project) and
   the combinatorial route (class-size-weighted within-class variance of the
   f_d(β)-weighted neighborhood profile) agree exactly, as the identity says.
3. Old hypothesis H-sample is re-verified — and upgraded: the original
   observation was about entrywise N estimation; what actually matters is
   that the *variance functional* of Theorem 3 concentrates this fast.

## Caveats

- Stratified sampling assumes class membership is enumerable (true here via
  full cost arrays; at n > 24 one would sample bitstrings and bin by cost,
  adding a stratification-error term not measured here).
- Errors quoted are for one-layer leakage at fixed angles; accumulated Σλ
  along schedules will average down sampling noise further.

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/008_sampled-leakage-predictor/run.jl
```

Seed 20260611. Output: `results.csv`. Smoke test: prefix `E1_SMOKE=1`.
