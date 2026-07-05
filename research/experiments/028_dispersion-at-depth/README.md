# 028 — Does the dispersion niche survive depth?

**Question.** The proxy's one won regime (E022–E026) is high weight
dispersion — measured at p=1, holding at p=3. But on unweighted graphs
every proxy variant collapses at p=10–20 ramps while transfer excels
(E019). The (high dispersion × depth) cell is unmeasured: does the niche
extend to depth, or do the argmax-drift losses of E019 overwhelm it?

**Answer.** *Pending — submitted 2026-07-05.*

## Method

Pareto(α=1.5) and lognormal(σ=2.5) weights (the two deep-dispersion
points), (ER(0.5), 3-regular) at n=16, 10 instances, p ∈ {10, 20} linear
ramps on per-instance mean-weight-scaled 8⁴ endpoint grids. Methods:
binned exact/sampled proxy (K=64, S=10), transfer×3 from weighted n=10
sources (index transfer ≡ rescaled transfer), pooled universal — against
true ramp ceilings on the same grid (regret ≥ 0, asserted).

## Reproduce

```bash
cd research/experiments/028_dispersion-at-depth && sbatch run.sb
E28_SMOKE=1 julia --project research/experiments/028_dispersion-at-depth/run.jl  # smoke
```
