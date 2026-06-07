# Systematic TriangleProxy vs NormalProxy Evaluation on Non-ER Graphs

**Date:** 2026-04-02
**Paper section:** approximation-ratios
**Tags:** triangle-proxy, normal-proxy, barabasi-albert, watts-strogatz, erdos-renyi, consistent-NP, sampling
**Status:** complete

## Motivation
Prior experiments showed fitted proxies fail when N/P are inconsistent. This
experiment tests both proxy types with *consistent* N/P pairing (each proxy
uses its own P(c')), and adds the new sampled homodist (SampN+EmpP) method.
This is a head-to-head comparison of all available proxy approaches.

## Setup
- **Script:** `julia/paper_figures/systematic_proxy_evaluation.jl`
- **Parameters:** n_fit=10, n_source=9, n_target=12, 50 fit instances, 20 eval instances
- **Depths:** p=1,2,3,5 (linear ramp for p>1)
- **Graph types:** ER(p=0.5), BA(m=2), WS(k=4,p=0.3)
- **Methods compared:**
  1. Transfer (median params from small source graphs)
  2. PaperProxy (analytical N+P, with p_eff for non-ER)
  3. Tri+TriP (fitted TriangleProxy N + TriangleProxy P)
  4. Norm+NormP (fitted NormalProxy N + NormalProxy P)
  5. EmpN+EmpP (exact empirical N + empirical P, multi-instance averaged)
  6. SampN+EmpP (sampled N with S=10 + empirical P)
- **Output files:**
  - Figure: `julia/paper_figures/output/systematic_proxy_evaluation.png`

## Fitted Proxy Parameters

| Graph | TriangleProxy (h,c,la,ra) | MSE | NormalProxy (cm,c1,c2) | MSE |
|-------|---------------------------|-----|------------------------|-----|
| ER | 5.19, 0.32, 0.38, 0.13 | 5.2e-8 | 11.4, 4.2, 5.1 | 5.1e-8 |
| BA | 6.00, 0.70, 0.18, 0.38 | 1.9e-7 | 8.8, 2.9, 3.1 | 1.4e-7 |
| WS | 5.64, 0.64, 0.19, 0.32 | 1.3e-7 | 10.3, 3.1, 3.5 | 9.4e-8 |

Both proxies achieve low MSE against empirical homodist, but this doesn't
predict QAOA performance (consistent with np-consistency finding).

## Key Findings

### 1. TriangleProxy degrades catastrophically at p>1
At p=1, Tri+TriP is competitive with PaperProxy (~0.78). But by p=3-5 it
drops to 0.65-0.69 — worse than any other method. The triangle shape is too
crude to sustain accurate multi-layer QAOA propagation. Errors compound
across layers.

### 2. NormalProxy is consistently the worst proxy
NormalProxy achieves 0.62-0.79 across all configurations. At p=1, it's already
poor (0.63-0.67 for non-ER). The covariance matrix parameterization doesn't
capture the asymmetric, multi-modal structure of N(c';d,c).

### 3. EmpN+EmpP and SampN+EmpP are the best non-analytical proxy methods
These two are effectively tied, confirming that S=10 sampling faithfully
reproduces exact homodist. At p=2-3 on non-ER, SampN+EmpP slightly
outperforms EmpN+EmpP (e.g., WS p=3: 0.819 vs 0.806), likely due to
beneficial noise regularization.

### 4. PaperProxy with p_eff surprisingly strong on non-ER at high depth
At p=3-5 on BA/WS, PaperProxy outperforms EmpN+EmpP despite being designed
for ER:
- BA p=5: PaperProxy 0.877 vs EmpN+EmpP 0.809
- WS p=5: PaperProxy 0.860 vs EmpN+EmpP 0.819

The analytical formula's regularization effect dominates the model mismatch
for non-ER graphs. The smoothing from the multinomial formula helps the
optimizer find better parameters, even though the formula is technically
wrong for BA/WS.

### 5. Transfer dominates all proxy methods
Transfer consistently achieves 0.05-0.10 higher approximation ratios than
the best proxy method at every depth and graph type:
- ER p=5: Transfer 0.902 vs EmpN+EmpP 0.820
- BA p=5: Transfer 0.926 vs PaperProxy 0.877
- WS p=5: Transfer 0.929 vs PaperProxy 0.860

## Significance

**For the paper:**
- TriangleProxy and NormalProxy are NOT competitive replacements for
  PaperProxy, even with consistent N/P pairing. They work at p=1 but
  degrade at higher depths.
- The alternative proxy shapes (Contribution 2 from our research plan) are
  a **negative result**: simpler functional forms cannot capture enough
  structure in N(c';d,c) for multi-layer propagation.
- EmpN+EmpP (and SampN+EmpP) provide the best non-analytical proxy, but
  even PaperProxy with p_eff beats them on non-ER at higher depths.
- The sampling estimator (Contribution 1) is validated and produces results
  equivalent to exact homodist.

**Implications for research direction:**
- The proxy approach itself (not just proxy shapes) may be fundamentally
  limited compared to transfer. The ~0.07 gap is persistent and doesn't
  close with better proxies or more data.
- The paper should frame TriangleProxy/NormalProxy as negative results that
  demonstrate the difficulty of approximating N(c';d,c) with simple shapes.
- The sampling estimator is the main positive methodological contribution.

## Next Steps Arising
- [ ] Investigate why PaperProxy+p_eff outperforms empirical homodist on non-ER
  at high depth — is it the smoothing/regularization, or is p_eff accidentally
  accurate for BA/WS cost structure?
- [ ] Consider proxy+refinement as a viable approach: use proxy (any type) to
  warmstart, then refine on target graphs. This was shown to help in
  hybrid-warmstart entry.
