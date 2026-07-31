# E025 — Does the recommended recipe compose? Sampled N + normalized objective

**Question**: Does the paper's recommended recipe — sampled N (S=10 per cost
class) *combined with* the normalized objective at depth — actually deliver
the regret its two separately-validated ingredients (E010: sampled N,
unnormalized; E020/E024: normalization, exact N) promise?

**Answer**: **Yes, and better than either ingredient alone: sampled N (S=10)
+ normalized objective is the best of all four combinations at BOTH depths —
it beats even exact N + normalized (p=1 mean regret 0.021 vs 0.047, 123/7
win/loss, sign p=8e-29; p=3 0.024 vs 0.044, 112/23, p=2e-15), uniformly
across all 7 families.**

## Method

Standard instance set: 7 families × n ∈ {12, 14} × 10 instances, fixed seeds
(`instance_seed(fi, n, inst)`, base seed 20260611). For each instance and each
depth (p=1 grid; p=3 linear-ramp grid), compute the true QAOA ceiling by
statevector simulation, then the regret of the grid point selected by each of
the four proxy objectives: {exact N, sampled N (S=10)} × {unnormalized,
normalized}. Sampling uses an independent RNG stream (seed + 777).

## Result

Mean (median) regret over 140 instances:

| combination      | p=1           | p=3           |
|------------------|---------------|---------------|
| exact N, raw     | 0.031 (0.031) | 0.080 (0.080) |
| exact N, norm    | 0.047 (0.044) | 0.044 (0.036) |
| sampled N, raw   | 0.030 (0.030) | 0.073 (0.073) |
| sampled N, norm  | **0.021 (0.017)** | **0.024 (0.021)** |

The exact-N column reproduces E020/E024 exactly (p=3: 0.080 -> 0.044; p=1
normalization hurts with exact N). The surprise: with *sampled* N the
normalized objective helps at p=1 too, and the combination beats exact N +
normalized at both depths, in every family (per-family means in
`results.csv`; sign tests above). Sampling noise in N appears to act as a
regularizer against the proxy's winner's curse (E022): the systematic
overprediction at the proxy's own argmax is what selection exploits, and
S=10 estimation noise plus normalization disrupts it.

## Caveats

- Grid ceilings, not continuous optima (E018 showed these are tight).
- S=10 per class is the single sampling budget tested; no S-sweep here.
- One sampling draw per instance (rng seed+777). The 140 instances make the
  *average* effect solid, but whether the gain comes from noise per se or
  from a systematic bias of the S=10 estimator (shrinkage toward uniform?)
  is NOT distinguished here — that is the E026 mechanism question. Do not
  put the "sampled beats exact" claim in the paper until E026 resolves it.

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/025_sampled-plus-normalized/run.jl   # ~15 min
E25_SMOKE=1 julia --project research/experiments/025_sampled-plus-normalized/run.jl              # quick check
```
