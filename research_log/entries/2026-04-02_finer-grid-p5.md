# Finer Grid at p=5: Proxy Landscape Degrades at Higher Depth

**Date:** 2026-04-02
**Paper section:** discussion
**Tags:** grid-resolution, high-depth, proxy-accuracy, overfitting
**Status:** complete

## Motivation
At p=5, Transfer beats SampN+EmpP. One hypothesis: the 10^4 ramp grid is too
coarse. Test finer grids (15^4, 20^4) to see if resolution is the bottleneck.

## Setup
- **Script:** `julia/paper_figures/finer_grid_p5.jl`
- **Parameters:** n=18, p=5, ER(p=0.5), S=20, gs=10/15/20
- **Grid sizes:** 10^4=10000, 15^4=50625, 20^4=160000

## Key Findings

### Finer grids make SampN+EmpP WORSE

| Method | Grid size | QAOA ratio | Proxy exp | Time |
|--------|-----------|-----------|-----------|------|
| Transfer | 10^4 | 0.908 | N/A | - |
| Transfer | 15^4 | **0.916** | N/A | - |
| Transfer | 20^4 | 0.910 | N/A | - |
| SampN+EmpP | 10^4 | **0.882** | 92.7 | 9.6s |
| SampN+EmpP | 15^4 | 0.875 | 97.2 | 30s |
| SampN+EmpP | 20^4 | 0.871 | 97.5 | 98s |

Higher proxy expectation value (97.5 > 92.7) corresponds to LOWER real QAOA
performance (0.871 < 0.882).

### Interpretation: proxy landscape overfitting

At p=5, the proxy landscape (under perfect homogeneity) diverges from the real
QAOA landscape. The finer grid finds parameters that are better *according to the
proxy* but worse *in reality*. This is classic overfitting: the proxy model is
insufficiently accurate at p=5, and optimizing harder against it exploits the
model's errors.

At p=1-3, the proxy landscape is a good approximation, so finer grids would help.
At p=5, the approximation breaks down.

Transfer doesn't have this problem because it optimizes on real QAOA (on small
instances). The finer grid helps Transfer (0.908 → 0.916), confirming that
resolution does matter for real QAOA optimization.

## Significance

**For the paper:**
- Explains why SampN+EmpP fails at p=5: the proxy landscape is inaccurate at
  higher depth, and grid search exploits the inaccuracy
- The coarser grid (10^4) acts as implicit regularization — it can't fully exploit
  the proxy's errors, resulting in better real QAOA performance
- This strongly motivates SampN+EmpP + Refinement: use the proxy for warmstart
  (finding the right neighborhood) then refine on real QAOA (finding the exact
  optimum). The proxy is still useful for warmstart even at p=5.
- Key message: "proxy-based methods should use coarse grid search + real QAOA
  refinement, not fine grid search on the proxy alone"

---

*Autonomous overnight run, 2026-04-02*
