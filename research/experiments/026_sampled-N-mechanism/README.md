# E026 — Why does sampled N beat exact N? Noise vs estimator bias

**Question**: Is E025's surprise (sampled N + normalized objective beating
*exact* N + normalized at both depths) driven by the sampling noise itself
interacting with argmax selection, or by a systematic property of the S=10
estimator?

**Answer**: **PENDING** (full run in progress).

## Method

Same instance seeds as E025, n=12 only (7 families × 10 instances). The
sampler (`sampled_homogeneous_distribution`) averages exact n(x;d,c) rows
over ≤S members per cost class, so it is unbiased for the class mean; the
two candidate mechanisms are separated by:

- **A. Replicates**: R=10 independent S=10 draws per instance → the
  per-replicate regret spread, plus the regret of the replicate-averaged
  N̄ (noise suppressed, any bias kept). If N̄ regresses to exact-N regret,
  noise is the active ingredient.
- **B. S-sweep**: S ∈ {3, 10, 30, 100, 300}, one draw each. If noise is
  the ingredient, the gain fades as S grows (sampled N → exact N).

Both depths (p=1 grid, p=3 linear-ramp grid), both objectives recorded.

## Result

PENDING.

## Caveats

- n=12 only (E025 established the effect at n ∈ {12, 14}).
- Grid ceilings, not continuous optima (tight per E018).

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/026_sampled-N-mechanism/run.jl   # ~20 min
E26_SMOKE=1 julia --project research/experiments/026_sampled-N-mechanism/run.jl              # quick check
```
