# E043 — Synthetic-noise control: generic noise vs sampling structure

**Question**: Does generic entrywise noise on the exact N reproduce the
sampled-N regret gain (E041/E042), or is the structure of real member-profile
heterogeneity essential?

**Answer**: **Depth dissociates the two: at p=1 iid multiplicative noise
reproduces and even beats the sampled gain (σ=0.05–0.1 gives mean regret
0.011–0.019 vs sampled 0.023, exact 0.042; U-shaped in σ), but at p=3 iid
noise only hurts (0.036–0.126, all σ, vs sampled 0.0225, exact 0.0397) —
at depth the gain requires the sampler's structured noise, not noise per se.**

## Method

Same 70 instances (7 families × n=12, E041 seeds). Variants: exact N; one
S=10 sampled draw; N′ = N · max(0, 1 + σz) with iid z ~ Normal(0,1),
σ ∈ {0.05, 0.1, 0.2, 0.4, 0.8}, 5 replicates each (350 draws per σ per
depth). Both depths, both objectives.

## Result

Normalized-objective mean regret (raw in parentheses):

| variant  | p=1             | p=3             |
|----------|-----------------|-----------------|
| exact    | 0.042 (0.030)   | 0.040 (0.075)   |
| samp S=10| 0.023 (0.031)   | 0.023 (0.075)   |
| iid σ=0.05| 0.019 (0.142)  | 0.036 (0.238)   |
| iid σ=0.1 | **0.011** (0.157) | 0.043 (0.254) |
| iid σ=0.2 | 0.017 (0.179)  | 0.053 (0.290)   |
| iid σ=0.4 | 0.043 (0.188)  | 0.074 (0.295)   |
| iid σ=0.8 | 0.089 (0.192)  | 0.126 (0.299)   |

Two structural facts: (1) iid noise destroys the *raw* objective at every σ
and depth (0.14–0.30) while sampling noise leaves it untouched — the
sampler's rows are genuine profiles with exact marginals (each row sums to
binom(n,d)), so its noise lives inside the physical constraint set; iid
noise leaves it. (2) The normalized objective converts small in-set noise
into a selection benefit at both depths (sampled), but out-of-set noise only
survives this conversion at p=1, where a single layer touches N once.

Practical reading: the recipe's cheap estimator is not "noisy but good
enough" — its noise is the *right kind*, and adding synthetic noise is not a
substitute at depth.

## Caveats

- n=12 only; multiplicative Gaussian is one noise model (clipped at 0).
- 5 replicates per σ (350 draws pooled per cell); single S=10 reference draw
  per instance (E042 has 10).

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/043_synthetic-noise-control/run.jl   # ~15 min
E27_SMOKE=1 julia --project research/experiments/043_synthetic-noise-control/run.jl              # quick check
```
