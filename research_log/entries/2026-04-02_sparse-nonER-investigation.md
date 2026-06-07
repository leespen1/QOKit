# Sparse Non-ER Investigation

**Date:** 2026-04-02
**Script:** `julia/paper_figures/sparse_nonER_investigation.jl`
**Status:** Complete

## Question

Does SampN+EmpP work better on sparse non-ER graphs (BA(m=1), WS(k=2)) where
the random baseline is lower and there's more headroom for proxy improvement?

## Setup

- BA: m_attach ∈ {1, 2, 3}, WS: k ∈ {2, 4, 6} (p_rewire=0.3)
- n=12-18, p=1,3
- 20 homodist instances, 10 eval instances, S=20 samples/cost
- Methods: Random, Transfer, PaperProxy(p_eff), SampN+EmpP

## Key Results

**SampN+EmpP consistently loses to Transfer on ALL non-ER configurations:**

| Config | p_eff | Random | Transfer | PaperProxy | SampN+EmpP | Gap |
|--------|-------|--------|----------|------------|------------|-----|
| BA(m=1) n=18 p=1 | 0.11 | 0.500 | 0.709 | 0.702 | 0.662 | -0.048 |
| BA(m=1) n=18 p=3 | 0.11 | 0.500 | 0.837 | 0.804 | 0.728 | -0.109 |
| WS(k=2) n=18 p=1 | 0.12 | 0.521 | 0.762 | 0.753 | 0.713 | -0.049 |
| WS(k=2) n=18 p=3 | 0.12 | 0.521 | 0.875 | 0.819 | 0.768 | -0.107 |
| BA(m=2) n=18 p=1 | 0.22 | 0.612 | 0.777 | 0.766 | 0.730 | -0.048 |
| WS(k=4) n=18 p=1 | 0.24 | 0.639 | 0.809 | 0.795 | 0.763 | -0.046 |

**Sparsity does NOT help the proxy.** The gap between SampN+EmpP and Transfer
is ~0.04-0.05 at p=1 and ~0.09-0.11 at p=3, regardless of graph density.

## Additional Findings

1. **PaperProxy(p_eff) works well at low p_eff (< 0.35)**: It tracks Transfer
   closely on sparse BA/WS (within 0.01-0.03), suggesting the analytical
   formula is a reasonable approximation even for non-ER when the graph is sparse.

2. **PaperProxy catastrophically fails at p_eff > 0.5**: WS(k=6) at n=12 has
   p_eff=0.545, and PaperProxy drops to 0.489 (below random baseline 0.693).

3. **BA(m=1) random baseline is exactly 0.50**: Tree graphs have mean cut =
   max edges / 2, and c_opt = all edges cut. This is the theoretical minimum.

4. **SampN+EmpP gap worsens with n**: On BA(m=1), the gap goes from -0.035
   (n=12) to -0.048 (n=18) at p=1. On WS(k=2), from -0.037 to -0.049.

## Why SampN+EmpP Fails on Non-ER

The proxy heuristic's advantage on ER comes from the high degree of homogeneity
in ER graphs — bitstrings with the same cost truly behave similarly. Non-ER
graphs have:
- **Structured topology**: BA has hubs, WS has clusters. Bitstrings with the
  same cost can behave very differently depending on which specific vertices
  are cut.
- **Low effective n**: Sparse graphs have few edges, so the homodist has fewer
  data points to average over, increasing noise.
- **Stronger locality**: The cost landscape of non-ER graphs has more local
  structure that the cost-class averaging in homodist destroys.

## Implications

- **SampN+EmpP is an ER-specific method.** It should only be recommended for
  ER graphs (or perhaps other random graph families with high homogeneity).
- **For non-ER graphs, Transfer remains the best practical method** at all
  depths and densities tested.
- **PaperProxy(p_eff) is a viable proxy for sparse non-ER** (p_eff < 0.35)
  but not for dense non-ER.
- **The paper should clearly scope its proxy claims**: the sampling-based
  homodist contribution is for ER graphs; non-ER requires different approaches.
