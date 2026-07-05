# 028 — Does the dispersion niche survive depth?

**Question.** The proxy's one won regime (E022–E026) is high weight
dispersion — measured at p=1, holding at p=3. But on unweighted graphs
every proxy variant collapses at p=10–20 ramps while transfer excels
(E019). The (high dispersion × depth) cell is unmeasured: does the niche
extend to depth, or do the argmax-drift losses of E019 overwhelm it?

**Answer. The niche narrows at depth but does not vanish — and the
(dispersion × depth) corner is hard for every method.** Under lognormal
σ=2.5 at p=10/20 the binned proxies (0.05–0.11) beat every single
transfer source (0.07–0.29; sources are again a lottery), though the
pooled universal schedule recovers on 3-regular (0.035) and the proxy
leads outright on dense ER at p=20 (0.059 vs 0.087 best source). Under
Pareto α=1.5 at depth everything degrades toward parity (proxy
0.11–0.20, best sources 0.03–0.13 with huge spread) — no method is good.
Array 11770562, all COMPLETED (~6 min/task).

Reading: depth compounds the E019 argmax-drift for proxies AND breaks
source-to-target alignment for transfer; under dispersion the proxy
remains the most *reliable* single-shot choice (lognormal), but the deep
heavy-tail corner is honestly open — plausibly needing finer schedule
grids or adaptive methods beyond this paper's scope.

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
