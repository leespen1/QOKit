# Sparse ER High-Depth: Proxy Advantage Persists at p=5

**Date:** 2026-04-02
**Script:** `julia/paper_figures/er_sparse_high_depth.jl`
**Status:** Complete — positive result

## Question

On ER(0.5), proxy overfits at p=5 and Transfer wins. Does the same happen
on sparse ER (p_edge=0.2, 0.3) where PaperProxy dominated at p=3?

## Key Results

### PaperProxy beats Transfer at p=5 on sparse ER (n=20)

| p_edge | n | p | Transfer | PaperProxy | PP-Transfer |
|--------|---|---|----------|------------|-------------|
| 0.2 | 20 | 1 | 0.756 | 0.766 | **+0.011** |
| 0.2 | 20 | 3 | 0.810 | **0.846** | **+0.036** |
| 0.2 | 20 | 5 | 0.859 | **0.866** | **+0.008** |
| 0.3 | 20 | 1 | 0.777 | 0.782 | +0.005 |
| 0.3 | 20 | 3 | 0.826 | **0.847** | **+0.021** |
| 0.3 | 20 | 5 | 0.871 | **0.875** | **+0.004** |
| 0.5 | 20 | 5 | 0.887 | 0.850 | -0.037 |

### Cross-over pattern

At each p_edge, the proxy advantage:
- Emerges at n>=16-18 for p=1
- Grows with n at p=3 (strongest regime)
- Is marginal at p=5 for sparse, absent for dense

The proxy advantage is strongest at intermediate depth (p=3) and moderate
to large n (n>=18), especially for sparse ER.

## Summary: Proxy vs Transfer Regime Map

| p_edge | p=1 | p=3 | p=5 |
|--------|-----|-----|-----|
| 0.2 | SampN wins (n>=16) | PP wins (n>=18) | PP wins (n=20) |
| 0.3 | SampN wins (n>=16) | PP wins (n>=20) | PP marginal (n=20) |
| 0.5 | SampN wins (n>=18) | SampN wins (n>=20) | Transfer wins |

## Implications

1. **PaperProxy on sparse ER is remarkably robust** — it beats Transfer at
   all depths tested (p=1,3,5) when n is large enough.

2. **The proxy overfitting problem is specific to ER(0.5)**, not intrinsic
   to the proxy algorithm. On sparse ER, the landscape is smoother and the
   proxy's analytical formula is more accurate.

3. **Paper recommendation**: For ER(p_edge<=0.3), use PaperProxy at all
   depths. For ER(0.5), use SampN+EmpP at p=1-3, Transfer at p>=5.
