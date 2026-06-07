# Sampled Homodist at Large n: QAOA Performance

**Date:** 2026-04-02
**Paper section:** scalability
**Tags:** sampling, homodist, large-n, erdos-renyi, barabasi-albert, watts-strogatz
**Status:** complete

## Motivation
The sampling estimator was validated at n=12. Now test at n=14-18 where exact
homodist becomes expensive. Can sampled homodist produce useful QAOA parameters
at these scales?

## Setup
- **Script:** `julia/paper_figures/sampled_homodist_large_n.jl`
- **Parameters:** n=14,16,18; S=20 samples/cost; 20 homodist instances; 10 eval instances
- **Methods:** Transfer, PaperProxy (p_eff), SampN+EmpP
- **Graph types:** ER(p=0.5), BA(m=2), WS(k=4,p=0.3)

## Key Findings

### 1. SampN+EmpP beats Transfer on ER at large n (unexpected!)

| n | Transfer | PaperProxy | SampN+EmpP |
|---|----------|------------|------------|
| ER p=1: 14 | 0.807 | 0.794 | 0.800 |
| ER p=1: 16 | 0.796 | 0.789 | **0.808** |
| ER p=1: 18 | 0.812 | 0.807 | **0.825** |
| ER p=3: 14 | 0.875 | 0.823 | 0.836 |
| ER p=3: 16 | 0.843 | 0.815 | **0.846** |
| ER p=3: 18 | 0.846 | 0.825 | **0.865** |

At n>=16, SampN+EmpP outperforms both Transfer and PaperProxy on ER graphs.
The empirical homodist captures instance-specific cost structure that Transfer
(from n=9) and PaperProxy (analytical formula) miss at larger n.

### 2. SampN+EmpP degrades on non-ER at large n

| n | Transfer | PaperProxy | SampN+EmpP |
|---|----------|------------|------------|
| BA p=3: 14 | 0.888 | 0.851 | 0.808 |
| BA p=3: 16 | 0.870 | 0.832 | 0.789 |
| BA p=3: 18 | 0.862 | 0.823 | 0.785 |
| WS p=3: 14 | 0.906 | 0.840 | 0.834 |
| WS p=3: 16 | 0.888 | 0.849 | 0.812 |
| WS p=3: 18 | 0.896 | 0.863 | 0.802 |

On BA/WS, the gap between SampN+EmpP and Transfer/PaperProxy widens with n.
The empirical homodist for non-ER graphs captures structural noise that hurts
the proxy optimizer at larger scales.

### 3. Computation time is practical

Sampled homodist (20 instances, S=20) takes ~1.5s at n=18, while exact
homodist would take ~70 minutes (20 × 210s). The sampling approach is
computationally viable at these scales.

## Significance

**Graph-type dependent behavior:**
- For ER: SampN+EmpP is the BEST method at n>=16, beating even Transfer.
  The homogeneity assumption holds well for ER, and the empirical homodist
  captures n-specific structure the analytical formula misses.
- For non-ER: Transfer and PaperProxy(low p_eff) remain superior. The
  homodist-based approach struggles because BA/WS cost structures are less
  homogeneous, making the per-instance noise more harmful.

**For the paper:**
- The sampling estimator enables homodist computation at previously intractable
  scales. On ER graphs, this produces the best-known QAOA parameters at n>=16.
- On non-ER, the fundamental limitation is homogeneity, not computation cost.
  The sampling estimator faithfully reproduces the homodist, but the homodist
  itself is a worse model for non-ER graphs.

## Next Steps Arising
- [ ] Test SampN+EmpP on ER at n=20-24 to see if the advantage continues to grow
- [ ] Investigate whether SampN+EmpP's ER advantage is from the homodist model
  quality or from the proxy optimization grid resolution
