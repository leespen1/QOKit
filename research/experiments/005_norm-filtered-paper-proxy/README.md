# 005 — Does the norm filter make the analytical PaperProxy safe everywhere?

**Question:** Across all 420 instances of experiment 004, does rejecting grid
points where the analytical proxy's state norm inflates (a model-error
certificate by Theorem 1) repair PaperProxy's catastrophic ER(0.5) failures
without harming families where it already worked?

**Answer: PENDING — full run queued.**

## Why this matters

Experiment 004 found PaperProxy's argmax on dense ER(0.5) lands on unphysical
norm-inflated peaks (predicted ⟨C⟩ > number of edges), while it beats the exact
compression on every other family. The exact compression is contractive
(‖φ‖ only decays — Theorem 1), so model-induced norm inflation is a free
artifact detector. The open question from the smoke run: mild norm drift is
ubiquitous on sparse graphs, where a *strict* filter (‖φ‖² ≤ 1.05) over-filters
and picks bad parameters — so we compare strict vs loose (‖φ‖² ≤ 2, gross
inflation only) thresholds. If the loose filter repairs ER(0.5) and leaves the
rest untouched, "analytical N + norm sanity check" becomes the paper's
recommended practical recipe.

## Method

Same instances and grids as experiment 004. For the PaperProxy arm only, AR at
the argmax of the proxy objective restricted to grid points with model-internal
norm² ≤ tol, for tol ∈ {1.05 (strict), 2.0 (loose)}, vs unrestricted (raw);
graceful fallback to the least-inflated points when the filter empties the grid
(flagged by frac_sane). `diagnose_paper_artifact.jl` in experiment 004 documents
the single-instance anatomy.

## Result

*(pending)*

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/005_norm-filtered-paper-proxy/run.jl
```

Seed 20260611 (same instance set as experiments 002/004). Output: `results.csv`.
Smoke test: prefix `E1_SMOKE=1`.
