# E042 — Why does sampled N beat exact N? Noise vs estimator bias

**Question**: Is E041's surprise (sampled N + normalized objective beating
*exact* N + normalized at both depths) driven by the sampling noise itself
interacting with argmax selection, or by a systematic property of the S=10
estimator?

**Answer**: **Noise is the active ingredient, and it acts only through the
normalized objective: the replicate-averaged N loses the gain (regret back to
exact-N level), the S-sweep is monotone with SMALLER S better (S=3 mean
regret 0.012/0.016 at p=1/3 vs exact-norm 0.042/0.040, converging to exact by
S=300), and noise does nothing for the raw objective at any S.**

## Method

Same instance seeds as E041, n=12 only (7 families × 10 instances). The
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

Mean regret over 70 instances (normalized objective):

| variant             | p=1    | p=3    |
|---------------------|--------|--------|
| exact N             | 0.0416 | 0.0397 |
| S=10 (700 draws)    | 0.0212 | 0.0236 |
| replicate-averaged N| 0.0360 | 0.0372 |
| S=3                 | 0.0120 | 0.0160 |
| S=30                | 0.0305 | 0.0309 |
| S=100               | 0.0359 | 0.0387 |
| S=300               | 0.0421 | 0.0397 |

Raw-objective regret is flat across all variants (~0.030 at p=1, ~0.075 at
p=3): noise helps *only* under normalization. Per-instance replicate SD of
regret is ~0.011 (median), so single draws vary, but the mean effect is
unambiguous (700 independent draws). Interpretation: subsampling noise
disrupts the selection effect (E038's winner's curse) that the normalized
proxy's systematic errors feed; less data gives better parameters. Whether
generic entrywise noise reproduces this, or the specific structure of real
member-profile heterogeneity is needed, is E043.

## Caveats

- n=12 only (E041 established the effect at n ∈ {12, 14}).
- Grid ceilings, not continuous optima (tight per E034).

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/042_sampled-N-mechanism/run.jl   # ~20 min
E26_SMOKE=1 julia --project research/experiments/042_sampled-N-mechanism/run.jl              # quick check
```
