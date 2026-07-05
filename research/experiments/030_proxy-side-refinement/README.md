# 030 — Free refinement: does polishing on the proxy's landscape transfer?

**Question.** E029 showed continuous refinement on the *true* landscape
erases the deep high-dispersion corner — at ~2000 true evaluations,
unaffordable where evaluations cost shots. The binned proxy's landscape
is continuous and free. If compass refinement on the *proxy's* predicted
objective transfers its gains to the true objective, the proxy earns a
new role: free continuous polish where evaluations are precious.

**Answer.** *Pending — submitted 2026-07-05.*

## Method

E028/E029's two worst cells (pareto15 × ER(0.5), logn2.5 × 3-regular;
n=16, p=20, same seeds, 5 instances). Start from the binned-sampled
proxy's grid choice; compass-refine ramp endpoints, then all 2p angles,
maximizing the *proxy's* predicted value (K=64, S=10); spend one true
evaluation per stage on the result. Compare true AR at: grid choice,
proxy-refined ramp4, proxy-refined full2p — reported against the 8⁴ grid
ceiling (E029's refined V* provides outside context: true-refinement
reaches ≈0.00 residual).

## Reproduce

```bash
cd research/experiments/030_proxy-side-refinement && sbatch run.sb
E30_SMOKE=1 julia --project research/experiments/030_proxy-side-refinement/run.jl  # smoke
```
