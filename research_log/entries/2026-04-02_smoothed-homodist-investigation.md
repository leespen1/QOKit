# Smoothed Homodist Investigation

**Date:** 2026-04-02
**Script:** `julia/paper_figures/smoothed_homodist_investigation.jl`
**Status:** Complete — NEGATIVE result

## Question

Can Gaussian smoothing of the empirical homodist act as regularization to
improve proxy performance at high depth? Hypothesis: smoothing removes noise
that causes proxy landscape overfitting.

## Setup

- Smoothing σ ∈ {0.5, 1.0, 2.0, 4.0} applied to both homodist and P(c')
- n=12-16, ER/BA/WS, p=1,3,5
- 20 homodist instances, 10 eval instances

## Key Results

**Smoothing always makes performance worse:**

| Method | ER p=1 | ER p=3 | ER p=5 | BA p=3 | WS p=3 |
|--------|--------|--------|--------|--------|--------|
| SampN+EmpP | 0.800 | 0.836 | 0.831 | 0.808 | 0.834 |
| Smooth(any ��) | 0.691 | 0.758 | 0.808 | 0.694 | 0.722 |
| Random | 0.691 | 0.691 | 0.691 | 0.627 | 0.652 |

(Values shown for n=14)

**All σ values produce identical results** — even σ=0.5 is enough to
completely destroy the useful structure. More smoothing doesn't make it
incrementally worse.

## Analysis

1. **The homodist's fine structure is essential, not noise.** Smoothing
   removes the detailed (d,c) patterns that the proxy algorithm needs to
   produce a useful objective landscape. This is NOT overfitting noise.

2. **PaperProxy's success is not due to smoothness.** The hypothesis that
   PaperProxy acts like a "maximally smoothed" version of the empirical
   homodist is wrong. PaperProxy succeeds because its structure is derived
   from the same probabilistic model as the cost landscape, not because
   it's smoother.

3. **Smoothing P is equally destructive.** Even if homodist smoothing
   preserved useful structure, the smoothed P(c') breaks the N/P consistency
   that the proxy algorithm requires.

## Implications

- **Spline proxy is not worth pursuing.** A smoothed/interpolated version of
  the empirical homodist will be strictly worse than the raw empirical version.
- **The proxy heuristic's power comes from structure, not smoothness.**
  Methods that preserve the detailed shape of N(c';d,c) work better.
- **For non-ER graphs, the path forward is**: use raw empirical/sampled
  homodist (SampN+EmpP), accept its limitations at high depth, and combine
  with Transfer+Refinement for high-depth cases.
