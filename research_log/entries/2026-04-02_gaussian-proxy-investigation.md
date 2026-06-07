# Gaussian Blob Proxy Investigation

**Date:** 2026-04-02
**Script:** `julia/paper_figures/gaussian_proxy_investigation.jl`
**Status:** Complete — NEGATIVE result

## Question

Can a Gaussian blob proxy (more expressive than Triangle, 7 parameters) produce
competitive QAOA parameters on non-ER graphs? Does a richer functional form
survive multi-layer propagation?

## Setup

- GaussianProxy: N(c';d,c) modeled as Gaussian in c for each (c',d) pair
  - 7 tunable parameters: center_target, sigma_base, sigma_scale, height_base,
    height_scale, height_power, center_bias
- Fitted via random search (200 random + 500 perturbation) to minimize
  normalized homodist MSE
- Tested with consistent GaussianP and empirical P
- n=14, ER/BA/WS, p=1,3,5

## Key Results

**GaussianProxy performs near or below random baseline on all configurations:**

| Method | ER p=1 | ER p=3 | BA p=1 | BA p=3 | WS p=1 | WS p=3 |
|--------|--------|--------|--------|--------|--------|--------|
| Random | 0.691 | 0.691 | 0.627 | 0.627 | 0.652 | 0.652 |
| Transfer | 0.807 | 0.875 | 0.792 | 0.888 | 0.818 | 0.906 |
| SampN+EmpP | 0.794 | 0.828 | 0.756 | 0.808 | 0.788 | 0.834 |
| GaussFit+GaussP | 0.621 | 0.690 | 0.583 | 0.617 | 0.623 | 0.631 |
| GaussFit+EmpP | 0.602 | 0.697 | 0.586 | 0.617 | 0.632 | 0.631 |

## Analysis

The Gaussian proxy achieves low homodist MSE (good shape fit) but produces
terrible QAOA parameters. This confirms a pattern seen with TriangleProxy:

**Homodist MSE ≠ QAOA performance.** A proxy can match the distribution shape
well but still produce a proxy objective landscape whose optima don't transfer
to real QAOA. The problem is that small shape errors get amplified through
the multi-layer proxy algorithm (repeated matrix-vector multiplications).

The consistent Gaussian P doesn't help either — GaussFit+GaussP and
GaussFit+EmpP perform similarly poorly.

## Why Parametric Proxies Fail

1. **Error amplification**: Each proxy layer multiplies the state by the
   homodist-derived matrix. Small per-element errors compound exponentially.
2. **Wrong optimization target**: Minimizing homodist MSE doesn't minimize
   QAOA landscape error. The landscape depends on the full product of matrices
   across all layers.
3. **PaperProxy succeeds because it's derived from the same model as the
   cost structure**, not because its shape is accurate. Its errors are
   *correlated* with the true structure in a way that preserves landscape topology.
4. **SampN+EmpP succeeds because it IS the empirical homodist** — no fitting
   step introduces systematic bias.

## Implications

- **Parametric proxy fitting is a dead end** for QAOA parameter setting.
  Both Triangle (4 params) and Gaussian (7 params) fail despite good homodist
  fit. More parameters won't help.
- **The path forward is either**: (a) use the empirical homodist directly
  (SampN+EmpP), (b) use a regularized version of the empirical homodist
  (spline smoothing), or (c) find an analytical formula for the specific
  graph family.
- **End-to-end optimization** (fit proxy to maximize real QAOA performance)
  might work but is expensive and defeats the purpose of the proxy heuristic.
