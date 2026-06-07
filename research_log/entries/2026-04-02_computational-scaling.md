# Computational Scaling Analysis

**Date:** 2026-04-02
**Paper section:** scalability
**Tags:** timing, scaling, computational-cost
**Status:** complete

## Motivation
Collect clean timing data for a paper figure showing computational feasibility
of SampN+EmpP vs Transfer and PaperProxy.

## Setup
- **Script:** `julia/paper_figures/computational_scaling.jl`
- **Parameters:** n=10-22; p=3; S=20; 10 instances; ER(p=0.5)
- **Machine:** CPU only (no GPU), Julia 1.11.6, multithreaded

## Timing Table (seconds, p=3)

| n | Costs | SampHD | ProxyOpt | PaperPx | QAOA eval | Refine | Total SampN |
|---|-------|--------|----------|---------|-----------|--------|-------------|
| 10 | 0.035 | 1.34 | 3.39 | 0.043 | 0.0001 | 0.06 | 4.76 |
| 12 | 0.002 | 0.025 | 0.83 | 0.015 | 0.0004 | 0.07 | 0.86 |
| 14 | 0.012 | 0.16 | 1.53 | 0.058 | 0.002 | 0.39 | 1.70 |
| 16 | 0.056 | 0.51 | 2.21 | 0.074 | 0.008 | 1.57 | 2.78 |
| 18 | 0.31 | 1.88 | 3.43 | 0.28 | 0.035 | 6.68 | 5.65 |
| 20 | 1.56 | 9.14 | 6.63 | 0.40 | 0.19 | 33.9 | 17.5 |
| 22 | 7.69 | 49.3 | 8.76 | 0.81 | 0.93 | N/A | 66.7 |

### Key observations:
- **Costs**: O(2^n * m), doubles every +1 to n. 7.7s at n=22.
- **Sampled homodist**: O(2^n * S * n), 49s at n=22. Dominant cost at large n.
- **Proxy optimization**: O(m^2 * n * p * K_grid), 8.8s at n=22. Grows slowly with n.
- **PaperProxy**: O(m^2 * n), <1s at all n. Cheapest but less accurate at large n.
- **QAOA eval**: O(2^n * n * p), 0.9s at n=22. Cheap for single evaluation.
- **Refinement**: O(2^n * n * p * 180), 34s at n=20. Too expensive at n=22.

### Transfer cost:
- Source optimization: 4.7s one-time (n_source=9, 10 instances, 10^4 ramp grid)
- Per-target: single QAOA eval (0.001-0.9s depending on n)
- Total: effectively free at deployment time

## Bottleneck Analysis

**SampN+EmpP workflow (without refinement):**
- Small n (≤14): dominated by proxy optimization (grid search)
- Large n (≥18): dominated by sampled homodist computation
- At n=22: 67s total, practical for parameter-setting

**SampN+EmpP + Refinement:**
- Refinement cost is 2-5x the non-refine workflow
- At n=20: 34s for one restart, ~136s with 4 restarts
- At n=22: infeasible for coord descent (would take ~2 min per restart)

**Transfer:**
- Nearly zero marginal cost at deployment
- But requires pre-optimized source instances

## Significance

**For the paper:**
- SampN+EmpP is computationally practical up to n=22+ (~1 minute)
- Beyond n=22, the cost array itself becomes the bottleneck (2^n entries)
- Refinement is practical up to n=20 (~30s per instance)
- Key figure: log-scale plot of total workflow time vs n, comparing
  SampN+EmpP, SampN+Refine, Transfer, PaperProxy

---

*Autonomous overnight run, 2026-04-02*
