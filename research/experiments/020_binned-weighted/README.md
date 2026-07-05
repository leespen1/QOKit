# 020 — Binned compression for weighted MaxCut

**Question.** With generic real edge weights every cost class is a
singleton and the integer-cost compression trivializes (the paper's
top referee risk). Theorems 2–3 hold for *any* fixed partition — so does
the calculus survive if the partition is K quantile bins of the weighted
cost, and is the binned proxy still a usable parameter setter?

**Answer. Yes — the calculus and the proxy transfer to weighted MaxCut
wholesale: quantile-binned leakage obeys λ² ≈ λ_struct² + O(1/K), and at
K ≈ 64–128 the binned proxy's regret matches the unweighted
exact-compression levels almost exactly.**

Jobs 11744099 (tasks 1, 2, 4) + 11745830 (task 3 rerun after two Julia-1.12
GC aborts; fixed by an allocation-free sweep, results bit-identical on the
smoke). Full tables: [analysis.txt](analysis.txt). Highlights:

- **Leakage:** at small angles λ falls like ~1/√K (λ²·K constant until the
  structural floor emerges, e.g. 3-regular n=14 γ=0.2: λ²·K =
  0.010→0.005 over K=4→128); at large angles the structural compression
  floor dominates and binning refines it only marginally (ER(0.5) n=16
  γ=1.0: 0.62→0.50). Exactly the predicted two-term decomposition — the
  binning term is a *variance* addition, hence 1/K in λ².
- **Parameter setting:** regret improves monotonically in K and saturates
  by K ≈ 64: ER(0.5) n=16 p=1 reaches 0.038 (unweighted analogue ~0.032),
  3-regular n=16 p=3 reaches 0.103 (unweighted 0.105). ~100 bins recover
  integer-cost behavior; even K=8 is already within 2× of saturation.
- The K-dimensional binned proxy costs O(K²n) per layer per schedule —
  independent of m — so the weighted extension is *cheaper* to run than
  the integer-cost proxy on dense graphs.

## Method

ER(0.5) and 3-regular at n ∈ {14, 16}, 10 instances each, i.i.d. U[0,1]
edge weights (fixed seeds). For each K ∈ {4, 8, 16, 32, 64, 128}:

- **Leakage:** one true weighted QAOA layer from |+⟩, then project onto the
  K equal-population bins; λ recorded at three representative (γ, β). The
  prediction: a structural compression term plus a binning term that
  shrinks with the within-bin cost spread (~1/K).
- **Parameter setting:** the binned proxy — bin-label homodist through the
  existing integer machinery, bin-mean phases cis(−γc̄_b/2), expectation
  with bin sizes and bin means — chooses schedules on the same grids as the
  true weighted ceiling (p=1: 40×40 with γ ∈ [0, 2π], since mean weight 1/2
  halves the cost scale; p=3 ramps: 8⁴ with γ ∈ [0.1, 3.2]). Regret vs the
  true ceiling, nonnegative by construction (asserted).

Output: `results_task<ID>.csv`, long format with `kind` ∈ {ceiling,
leak:γ,β, binned_proxy}, columns `K, p, value, ceil, regret`.

## Caveats

- Quantile (equal-population) bins are one choice; equal-width bins would
  weight the tails differently. Quantile bins make every bin's statistics
  well-sampled, which is what the homogeneity average needs.
- The binned proxy's expectation uses bin-mean costs — an additional
  homogeneity-style approximation entering only through the argmax.
- U[0,1] weights only; no negative or heavy-tailed weights tested.

## Reproduce

```bash
cd research/experiments/020_binned-weighted && sbatch run.sb
E20_SMOKE=1 julia --project research/experiments/020_binned-weighted/run.jl  # smoke
```
