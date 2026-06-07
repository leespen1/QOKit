# Fitted TriangleProxy Evaluation on BA/WS Graphs

**Date:** 2026-04-01
**Paper section:** approximation-ratios
**Tags:** triangle-proxy, fitting, barabasi-albert, watts-strogatz, erdos-renyi, parameter-setting, homodist
**Status:** complete

## Motivation
PaperProxy(p_eff) was shown (2026-03-17) to fall behind parameter transfer at
depth > 1 for BA/WS graphs. TriangleProxy fitted to empirical N(c';d,c) data
could potentially do better by capturing structure-dependent corrections.

## Setup
- **Script:** `julia/paper_figures/fitted_triangle_investigation.jl`
- **Fitting:** n=10, 50 instances per graph type, averaged empirical homodist
- **Fitting method:** Grid search (1728 configs) + local refinement (smart random search, 3000 iters)
- **Proxy type:** IntuitiveTriangleProxy with 4 parameters (height_adjustment, center_adjustment, left_angle, right_angle)
- **Evaluation:** n_source=9, n_target=12, 20 instances, p=1,2,3
- **Methods compared:** Transfer (coord descent p=1, linear ramp p>1), PaperProxy(p_eff), FittedTriangleProxy
- **Graph params:** BA(m_attach=2), WS(k=4, p_rewire=0.3), ER(p_edge=0.5)
- **Output files:** `julia/paper_figures/output/fitted_triangle_comparison.png`, `fitted_triangle_heatmaps.png`

## Fitted Parameters
| Graph | height | center | left_angle | right_angle | MSE |
|-------|--------|--------|------------|-------------|-----|
| ER | 5.622 | 0.400 | 0.289 | 0.147 | 5.17e-08 |
| BA | 2.972 | 0.514 | 0.269 | 0.251 | 1.94e-07 |
| WS | 4.621 | 0.446 | 0.303 | 0.199 | 1.28e-07 |

## Key Findings

### 1. Homodist fit quality: FittedTriangle >> PaperProxy
FittedTriangleProxy achieves dramatically lower normalized MSE against empirical
homodist than PaperProxy(p_eff):
- BA: 10.3x lower MSE
- WS: 4.8x lower MSE
- ER: 0.9x (PaperProxy slightly better, as expected)

### 2. QAOA performance: FittedTriangle << PaperProxy (negative result)
Despite much better homodist fit, FittedTriangleProxy produces **worse** QAOA
approximation ratios, and the gap widens dramatically with depth:

| Graph | p | Transfer | PaperProxy | FittedTriangle |
|-------|---|----------|------------|----------------|
| BA | 1 | 0.793±0.034 | 0.784±0.033 | 0.780±0.033 |
| BA | 2 | 0.857±0.030 | 0.767±0.031 | 0.766±0.027 |
| BA | 3 | 0.891±0.027 | 0.847±0.031 | **0.648±0.020** |
| WS | 1 | 0.802±0.027 | 0.792±0.025 | 0.787±0.025 |
| WS | 2 | 0.862±0.024 | 0.797±0.022 | 0.732±0.021 |
| WS | 3 | 0.894±0.024 | 0.828±0.023 | **0.647±0.019** |

Win counts (FittedTriangle vs PaperProxy): 0/20 for BA p=1, 9/20 for BA p=2,
0/20 for BA p=3, 0/20 for WS all depths.

### 3. Pearson correlation: PaperProxy > FittedTriangle
Despite lower MSE, FittedTriangleProxy has lower Pearson correlation with
empirical distributions:
- BA: PaperProxy 0.9905 vs FittedTriangle 0.9526
- WS: PaperProxy 0.9914 vs FittedTriangle 0.9522

## Interpretation

**Internal consistency trumps fit quality.** This result extends the P(c')
investigation finding to the full proxy system:

1. **PaperProxy's consistency advantage**: Both N(c';d,c) and P(c') are derived
   from the same probabilistic model (binomial/multinomial over independent edges).
   The amplitudes Q computed from N are calibrated to P's assumptions. Even when
   applied to non-ER graphs with p_eff, this internal consistency is preserved.

2. **FittedTriangle's consistency problem**: The fitted N(c';d,c) matches
   empirical data well, but TriangleProxy's P(c') is a crude triangular
   distribution — not matched to the fitted N. This P/N mismatch causes the
   proxy's objective function to diverge from the true QAOA landscape, especially
   at higher depths where errors compound.

3. **Depth compounding**: At p=1, the gap is small (~0.005). At p=3, FittedTriangle
   degrades catastrophically to ~0.65 (near random). This suggests small P/N
   inconsistencies amplify multiplicatively across QAOA layers.

4. **Generalization across n**: The proxy was fitted on n=10 but evaluated on n=12.
   PaperProxy's analytical formula generalizes naturally across n because it's
   parametric. TriangleProxy's fitted parameters may not transfer well.

## Implications for Research Directions

- **Direction 1 (better fitted proxy shapes)**: Fitting N alone is insufficient.
  Any fitted proxy must also have a consistent P(c'). Future work should either:
  (a) derive P(c') analytically from the fitted N(c';d,c), or
  (b) fit N and P jointly to ensure consistency, or
  (c) use PaperProxy's binomial P(c') paired with fitted N.

- **Direction 2 (non-ER graphs)**: PaperProxy(p_eff) remains the best proxy
  approach even for non-ER graphs, despite being derived for ER. The effective
  edge probability adaptation is surprisingly robust.

- **Direction 5 (multi-instance averaged homodist)**: Same consistency issue would
  apply. An averaged empirical N paired with the wrong P will underperform.

## Next Steps Arising
- [ ] Test whether pairing fitted N(c';d,c) with PaperProxy's binomial P(c') resolves the consistency issue
- [ ] Investigate whether PaperProxy(p_eff) can be improved by adjusting p_eff beyond m/C(n,2)
- [ ] Consider deriving P(c') from fitted N via normalization: P(c') ~ sum_d N(c';d,c) / 2^n

---

*Autonomous overnight run, 2026-04-01*
