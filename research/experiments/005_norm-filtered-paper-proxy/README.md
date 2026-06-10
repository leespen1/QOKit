# 005 — Does the norm filter make the analytical PaperProxy safe everywhere?

**Question:** Across all 420 instances of experiment 004, does rejecting grid
points where the analytical proxy's state norm inflates (a model-error
certificate by Theorem 1) repair PaperProxy's catastrophic ER(0.5) failures
without harming families where it already worked?

**Answer: No — the norm filter helps only on dense ER(0.5) (strict tolerance
halves its regret, 0.159 → 0.077 at p=1) and is catastrophic on every other
family (regret 0.15–0.31 vs 0.01–0.08 raw), because the analytical model's
norm calibration is grossly wrong on sparse graphs even where its argmax is
nearly perfect. NEGATIVE RESULT for the universal-recipe idea; the argmax
location is far more robust to model error than the objective values.**

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

Mean regret (ceiling − AR), n=12+14 pooled, 60 instances/family:

| family        | raw p1 | strict p1 | loose p1 | frac sane p1 | raw p3 | strict p3 | loose p3 |
|---------------|--------|-----------|----------|--------------|--------|-----------|----------|
| ER(0.5)       | 0.159  | **0.077** | 0.108    | 0.63         | 0.128  | 0.181     | 0.190    |
| ER(0.25)      | 0.011  | 0.242     | 0.238    | 0.05         | 0.044  | 0.302     | 0.302    |
| BA(k=2)       | 0.012  | 0.243     | 0.216    | 0.08         | 0.035  | 0.281     | 0.281    |
| BA(k=4)       | 0.012  | 0.071     | 0.028    | 0.71         | 0.076  | 0.147     | 0.133    |
| WS(k=4;b=0.1) | 0.015  | 0.160     | 0.173    | 0.19         | 0.056  | 0.260     | 0.260    |
| WS(k=4;b=0.5) | 0.010  | 0.154     | 0.163    | 0.19         | 0.051  | 0.273     | 0.273    |
| 3-regular     | 0.010  | 0.226     | 0.237    | 0.03         | 0.043  | 0.311     | 0.311    |

- On sparse families, *no* grid point passes even the loose tolerance at p=3
  (frac_sane = 0) — the analytical norm calibration is broken everywhere there —
  yet the raw argmax achieves regret 0.01–0.08. Norm inflation does NOT
  co-locate with bad parameters; filtering by it actively selects bad ones.
- On ER(0.5) the strict filter genuinely helps (worst case 0.40 → 0.29) but
  remains ~2.5× worse than the exact compression (0.030, exp 004).
- **Refined picture:** model error destroys the analytical proxy's absolute
  calibration (values, norms) long before it moves the argmax. The dense-ER
  failure mode is a *competing spurious peak* overtaking the true one, not a
  shifted peak. The right guard must therefore reject only certainly-impossible
  predictions rather than all decalibration → experiment 006 tests the
  physicality filter (predicted ⟨C⟩ ≤ m).

## Caveats

- Tolerances tested were only {1.05, 2.0}; intermediate values cannot rescue
  the sparse families given frac_sane ≈ 0 at 2.0.
- The model-internal norm uses the binomial P; an empirical-P norm could
  calibrate differently (untested).

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/005_norm-filtered-paper-proxy/run.jl
```

Seed 20260611 (same instance set as experiments 002/004). Output: `results.csv`.
Smoke test: prefix `E1_SMOKE=1`.
