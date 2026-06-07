# BA/WS Parameter Sweep: Proxy and Transfer Robustness

**Date:** 2026-04-01
**Paper section:** approximation-ratios
**Tags:** barabasi-albert, watts-strogatz, parameter-sweep, robustness, paper-proxy, parameter-transfer, density
**Status:** complete

## Motivation
The smallworld investigation (2026-03-17) tested only BA(m=2) and WS(k=4,p=0.3).
This sweep tests whether results hold across the full parameter space.

## Setup
- **Script:** `julia/paper_figures/parameter_sweep_investigation.jl`
- **BA sweep:** m_attach ∈ {1, 2, 3, 4}
- **WS sweep:** (k, p_rewire) ∈ {(2, 0.1), (4, 0.3), (6, 0.5)}
- **Evaluation:** n_source=9, n_target=12, 20 instances, p=1,2,3
- **Methods:** Transfer (coord descent p=1, linear ramp p>1) vs PaperProxy(p_eff)
- **Output:** `julia/paper_figures/output/parameter_sweep_ba.png`, `parameter_sweep_ws.png`

## Edge Densities (n=12)
| Config | Avg edges | p_eff |
|--------|-----------|-------|
| BA m=1 | 11 | 0.167 |
| BA m=2 | 21 | 0.318 |
| BA m=3 | 30 | 0.455 |
| BA m=4 | 38 | 0.576 |
| WS k=2,p=0.1 | 12 | 0.182 |
| WS k=4,p=0.3 | 24 | 0.364 |
| WS k=6,p=0.5 | 36 | 0.545 |

## Key Findings

### 1. PaperProxy catastrophically fails for dense non-ER graphs
For BA m=4 (p_eff=0.576) and WS k=6,p=0.5 (p_eff=0.545), PaperProxy drops to
near-random performance:

| Config | p=1 Transfer | p=1 PaperProxy | Gap |
|--------|-------------|----------------|-----|
| BA m=4 | 0.836 | **0.478** | +0.358 |
| WS k=6,p=0.5 | 0.816 | **0.493** | +0.323 |

PaperProxy produces parameters that are **worse than random** for these dense
non-ER graphs. The gap persists at all depths (0.30 at p=3 for BA m=4).

### 2. Transfer is robust across all configurations
Transfer consistently achieves good approximation ratios regardless of graph
family or density:
- p=1: 0.72-0.84 across all configs
- p=3: 0.85-0.91 across all configs
- Standard deviation consistently 0.01-0.04

### 3. For sparse/moderate graphs, PaperProxy works reasonably
At lower densities (p_eff < 0.5), PaperProxy performs within 0.01-0.10 of
Transfer:
- BA m=1,2,3: Transfer-Proxy gap = 0.008-0.071 depending on depth
- WS k=2,4: Transfer-Proxy gap = 0.009-0.072 depending on depth

### 4. Density is the critical factor
The proxy failure correlates with graph density, not graph family. Both BA and
WS show catastrophic failure at similar p_eff (~0.55). Below p_eff~0.45, the
proxy works adequately.

## Interpretation

PaperProxy's analytical formula (Eq. 12-16) assumes edges are independent with
probability p_eff. For **dense** non-ER graphs, this assumption fails badly:
- Dense BA graphs have strong hub structure (degree distribution is power-law)
- Dense WS graphs have lattice-based correlations not captured by i.i.d. model
- The multinomial sum for N(c';d,c) produces a distribution that poorly matches
  the actual cost-distance structure, leading to a proxy landscape whose optima
  are far from the true QAOA optima

At lower densities, the central limit theorem effect makes all graph families
look more similar (the distribution approaches Gaussian), so the ER approximation
works better.

## Implications for the Paper

1. **PaperProxy(p_eff) cannot be recommended for dense non-ER graphs** (p_eff > 0.5).
   This is a significant limitation that should be documented.

2. **Parameter transfer is the more robust approach for non-ER graphs**, especially
   when graph density varies or is unknown.

3. The density threshold (~p_eff=0.5) gives practitioners a clear guideline:
   use PaperProxy for sparse graphs, transfer for dense ones.

4. This motivates developing proxy methods that account for degree heterogeneity
   (e.g., incorporating degree sequence information into the proxy formula).

## Next Steps Arising
- [ ] Investigate what happens at the boundary (p_eff ≈ 0.45-0.55) with finer resolution
- [ ] Test whether the density threshold varies with n (check at n=14,16,18)
- [ ] Consider degree-aware proxy formulations for dense non-ER graphs

---

*Autonomous overnight run, 2026-04-01*
