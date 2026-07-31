# E045 — Selection-aware sharpening of the regret certificate

**Question**: Can the E036 certificate's ~10x slack be reduced by exploiting
the E038 anatomy — specifically the proxy-computable margin term, which can
be subtracted with no additional inequality?

**Answer**: **The margin-subtracted certificate is rigorous, never worse,
and holds on 140/140 instance-depth pairs, but it closes little of the gap
(median tightness 10.2→9.5 at p=1, 11.3→10.6 at p=3; no row improves 2x):
the slack lives in the pointwise Cauchy–Schwarz step, not in ignoring the
selection structure.** A useful negative with one positive corollary: the
practitioner-uniform variant (no oracle θ*) is also valid on every row at
finite cost (mean bound ~1.0 vs the paper's oracle form ~0.47).

## Method

Derivation (no new inequality beyond E036's):
regret = [e(θ*) − e(θ̂)] − margin ≤ ε(θ*) + ε(θ̂) − margin, with
margin = F̂(θ̂) − F̂(θ*) ≥ 0 computed by the proxy itself. Variants: the
practitioner-uniform bound ε(θ̂) + max_θ[ε(θ) − margin(θ)] (removes θ*),
and the non-certified one-sided diagnostic ε(θ̂) − margin (valid iff
e(θ*) ≤ 0, which E038 observed at 277/280 points but nothing proxy-side
certifies). 7 families × n ∈ {12,14} × 5 instances, p=1 grid and p=3 ramps;
asserts ε ≥ |e| at every schedule of every landscape (pointwise machine
check of the ε lemma), the exact identity, margin ≥ 0, and validity of all
certified bounds on every row.

## Result

| quantity (median over 70 rows) | p=1 | p=3 |
|---|---|---|
| paper bound / regret | 10.18 | 11.33 |
| sharpened bound / regret | 9.49 | 10.57 |
| mean margin | 0.030 | 0.033 |
| mean paper bound | 0.478 | 0.457 |
| one-sided diagnostic valid | 70/70 | 70/70 (sign ok 69/70) |

The margin (~0.03) is an order of magnitude smaller than the ε sum (~0.46),
so subtracting it barely moves the ratio. The error-vector cosine between
θ* and θ̂ (median −0.43 at p=1, −0.05 at p=3) shows the two error vectors
are not strongly aligned, so a two-point correlated-error refinement has
limited headroom as well. Conclusion consistent with E038/E039: sharp
regret theory needs to model the signed selection effect, and no
Cauchy–Schwarz-style pointwise bound will get there by bookkeeping alone.

## Caveats

- 5 instances per family (70 rows per depth); single grid resolutions.
- The one-sided diagnostic's validity relies on overprediction at θ*,
  which is empirical (E038), not certified.

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/045_selection-aware-certificate/run.jl  # ~15 min
E45_SMOKE=1 julia --project research/experiments/045_selection-aware-certificate/run.jl             # quick check
```
