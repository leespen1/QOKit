# 021 — Does transfer's dominance survive random weights?

**Question.** Transfer crushed every proxy variant on unweighted families
(E018/E019) because instances within a family share near-identical optimal
angles. Random U[0,1] edge weights individualize instances — a weighted
source's angles need not fit a weighted target — while the binned proxy
(E020) reads each instance's own cost structure. Is the heterogeneous
weighted setting the proxy's first genuine win?

**Answer. No crossover yet — weights erode transfer (3–10×) but it still
wins every cell; only the zero-knowledge universal schedule collapses to
proxy level (dense ER, p=3: 0.069–0.085).**

Array 11748691, all COMPLETED (0:36–2:35 per task). Cell means:
best single-source transfer 0.003–0.022 vs binned proxies 0.033–0.111;
binned_sampled ≈ binned_exact (sampling regularizes again, echoing E010).
Direction confirmed: i.i.d. U[0,1] weights individualize instances enough
to hurt concentration-based methods — universal loses its unweighted
near-perfection entirely on dense ER at depth — but same-family weighted
sources retain distribution-level concentration and keep winning. The
crossover to per-instance surrogates, if it exists, needs stronger
heterogeneity (heavy tails, mixed scales, structured weights).

## Method

E020's weighted targets exactly (same seeds → same ceilings): (ER(0.5),
3-regular) × n ∈ {14, 16}, 10 instances, U[0,1] weights; p=1 grid
(γ ∈ [0, 2π]) and p=3 ramps. Methods, all choosing from the ceiling's own
grid: `binned_exact` (K=64 quantile bins, exact bin-label homodist),
`binned_sampled` (K=64, N from S=10 samples per bin — the O(S·K·2ⁿ)
practical variant), `transfer_1..3` (angles optimized on the true landscape
of an n=10 same-family instance *with its own random weights*), `universal`
(pooled argmax over the 6 weighted sources).

## Caveats

- One weight law (U[0,1]); mean-1/2 weights keep all instances at the same
  overall scale, which *helps* transfer — heavier-tailed or per-instance
  rescaled weights would stress it further.
- Transfer sources are weighted instances of the same size-10 family; a
  practitioner could also rescale angles by total weight, not tested here.

## Reproduce

```bash
cd research/experiments/021_weighted-transfer && sbatch run.sb
E21_SMOKE=1 julia --project research/experiments/021_weighted-transfer/run.jl  # smoke
```
