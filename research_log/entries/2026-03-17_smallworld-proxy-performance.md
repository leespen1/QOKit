# Small-World Graph Investigation: QAOA Proxy Performance on BA and WS Graphs

**Date:** 2026-03-17
**Paper section:** parameter-transfer
**Tags:** barabasi-albert, watts-strogatz, erdos-renyi, paper-proxy, parameter-transfer, depth-scaling, homogeneity
**Status:** complete

## Motivation
Investigate whether the QAOA proxy heuristic from the Parameter-Setting paper
extends to Barabasi-Albert (BA) and Watts-Strogatz (WS) graphs. The PaperProxy
formula is derived for ER graphs; we test it on non-ER types using an effective
edge probability p_eff = m / C(n,2).

## Setup
- **Script:** `julia/paper_figures/smallworld_investigation.jl`
- **Parameters:** n=10 (distributions), n_source=9, n_target=12, instances=20, depths=1-3, grid=10^4 points
- **Output files:**
  - Figures: `julia/paper_figures/output/smallworld_*.png`

## Key Findings
1. PaperProxy with p_eff achieves Pearson correlations >= 0.98 on BA/WS, but
   normalized MSE is 10x higher than ER (shape matches, scale doesn't).
2. Parameter transfer works well for all graph types (~0.80 at p=1 to ~0.89 at p=3).
3. PaperProxy falls behind transfer significantly at p>=2 (gap 0.04-0.09), even
   with identical linear ramp grid search. Gap is similar across graph types.
4. Cross-type parameter transfer works within 0.005 penalty at p=1.
5. BA/WS are 2-3x more homogeneous than ER (lower coefficient of variation).

## Significance
- The proxy's approximation degrades with depth due to compounding homogeneity
  error, not due to non-ER structure specifically (gap is similar for ER too).
- Cross-type transferability suggests QAOA landscape depends primarily on problem
  scale (n,m), not graph structure — important for practical parameter setting.
- Higher homogeneity of BA/WS means fitted proxies could work *better* on these
  types than on ER, motivating the TriangleProxy fitting direction.

## Next Steps Arising
- [x] Fit TriangleProxy/NormalProxy to non-ER empirical distributions (moved to next_steps P0)
- [x] Test robustness across BA/WS parameter space (moved to next_steps P0)
- [x] Scale to larger n (moved to next_steps P0)

---

*Original report: `ClaudeResearchTasks/Task1/smallworld_investigation.md`*
