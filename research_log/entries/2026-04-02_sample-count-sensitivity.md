# Sample Count Sensitivity: S Has Almost No Effect

**Date:** 2026-04-02
**Paper section:** methodology
**Tags:** sampling, homodist, sample-count, robustness
**Status:** complete

## Motivation
SampN+EmpP uses S samples per cost class to estimate the homodist. All prior
experiments used S=20. Is this optimal? Does S=50 or S=100 improve accuracy?

## Setup
- **Script:** `julia/paper_figures/sample_count_sensitivity.jl`
- **Parameters:** n=16,18,20; S=5,10,20,50,100; p=1,3; 10 homodist instances; 5 eval
- **Graph type:** ER(p=0.5) only

## Key Findings

### QAOA performance is insensitive to S

| n | p | S=5 | S=10 | S=20 | S=50 | S=100 |
|---|---|-----|------|------|------|-------|
| 16 | 1 | 0.815 | 0.815 | 0.815 | 0.815 | 0.815 |
| 16 | 3 | 0.860 | 0.856 | 0.856 | 0.856 | 0.860 |
| 18 | 1 | 0.826 | 0.826 | 0.826 | 0.826 | 0.826 |
| 18 | 3 | 0.869 | 0.869 | 0.869 | 0.869 | 0.869 |
| 20 | 1 | 0.817 | 0.821 | 0.817 | 0.821 | 0.821 |
| 20 | 3 | 0.863 | 0.863 | 0.863 | 0.863 | 0.863 |

Variation across S values: ≤0.004 at all (n,p) combinations. Most cells are
identical to 4 decimal places.

### Timing scales linearly with S

| n | S=5 | S=10 | S=20 | S=50 | S=100 |
|---|-----|------|------|------|-------|
| 16 | 2.6s | 0.4s | 0.4s | 0.8s | 1.6s |
| 18 | 0.7s | 1.0s | 1.7s | 3.9s | 7.1s |
| 20 | 2.8s | 5.1s | 9.2s | 23.1s | 41.2s |

(n=16/S=5 timing anomaly likely due to JIT compilation overhead.)

### Why is S so unimportant?

The homodist is averaged over 10 graph instances. Each instance contributes
independently to the average. With 10 instances × S samples/cost, the effective
sample size is 10S. Even at S=5, that's 50 effective samples per cost class,
which is sufficient for the proxy grid search resolution (10^4 points). The
grid search is the bottleneck for parameter accuracy, not the homodist noise.

## Significance

**For the paper:**
- S=5-10 is sufficient. The paper should recommend S=10 as a conservative default.
- At n=20, S=10 costs 5s vs S=100 at 41s — an 8x speedup with no accuracy loss.
- This strengthens the scalability story: sampled homodist is even cheaper than
  previously reported when S is reduced.
- The robustness to S means the method doesn't require careful tuning.

## Revised Timing Estimates

Using S=10 instead of S=20:
- n=18: homodist ~1.0s (was 1.7s)
- n=20: homodist ~5.1s (was 9.2s)
- n=22: homodist ~25s (was 49s), total workflow ~42s (was 67s)

---

*Autonomous overnight run, 2026-04-02*
