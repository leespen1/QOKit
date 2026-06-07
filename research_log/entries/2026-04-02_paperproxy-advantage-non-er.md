# Why PaperProxy+p_eff Outperforms Empirical Homodist on Non-ER Graphs

**Date:** 2026-04-02
**Paper section:** discussion
**Tags:** paper-proxy, effective-edge-probability, regularization, barabasi-albert, watts-strogatz
**Status:** complete

## Motivation
The systematic proxy evaluation showed PaperProxy with effective edge probability
outperforms EmpN+EmpP on BA/WS at p>=3, despite the analytical formula being
derived for ER only. This investigation tests why.

## Setup
- **Script:** `julia/paper_figures/paperproxy_advantage_investigation.jl`
- **Tests:**
  1. Vary p_eff for PaperProxy on non-ER (0.1 to 0.8)
  2. Vary number of homodist instances for EmpN+EmpP (5 to 100)
  3. PaperProxy N + Empirical P vs PaperProxy N + Binomial P
  4. Landscape heatmap comparison
- **Output:** `julia/paper_figures/output/paperproxy_advantage_landscapes.png`

## Key Findings

### 1. Low p_eff is better than natural p_eff for PaperProxy on non-ER

| p_eff | BA p=1 | BA p=3 | WS p=1 | WS p=3 |
|-------|--------|--------|--------|--------|
| 0.1 | 0.764 | **0.855** | 0.769 | **0.863** |
| 0.2 | 0.774 | **0.855** | 0.780 | **0.863** |
| 0.3 | 0.782 | 0.852 | 0.789 | 0.841 |
| 0.4 | **0.787** | 0.828 | **0.796** | 0.825 |
| 0.5 | 0.769 | 0.799 | 0.778 | 0.809 |
| 0.6 | 0.407 | 0.695 | 0.410 | 0.649 |
| natural | ~0.32 | ~0.32 | ~0.36 | ~0.36 |

At p=3, the optimal p_eff is 0.1-0.2, well BELOW the natural density (~0.32 for
BA, ~0.36 for WS). PaperProxy with low p_eff creates a very smooth, broad
N(c';d,c) that acts as a strong regularizer. The specific graph structure doesn't
matter — it's the smoothness that helps.

The catastrophic collapse at p_eff>=0.6 (confirmed from parameter-sweep results)
is a phase transition in the multinomial formula.

### 2. More homodist instances don't help EmpN+EmpP

| N_inst | BA p=1 | BA p=3 | WS p=1 | WS p=3 |
|--------|--------|--------|--------|--------|
| 5 | 0.760 | 0.809 | 0.764 | 0.800 |
| 10 | 0.760 | 0.809 | 0.769 | 0.806 |
| 20 | 0.760 | 0.809 | 0.769 | 0.806 |
| 50 | 0.760 | 0.809 | 0.769 | 0.806 |
| 100 | 0.760 | 0.809 | 0.769 | 0.806 |

The empirical homodist average saturates at ~10 instances. The gap with
PaperProxy (0.85 vs 0.81 at p=3) is NOT due to insufficient averaging.
The empirical N(c';d,c) is noisier/more structured than PaperProxy's smooth
formula, and this extra structure hurts the optimizer.

### 3. P(c') source doesn't matter for PaperProxy

PaperN+BinomialP vs PaperN+EmpiricalP are essentially identical (within 0.003).
The advantage is entirely in the N(c';d,c) computation, not in P(c'). This
is surprising given the N/P consistency principle — but it makes sense because
PaperProxy's N is so smooth that it works with any reasonable P.

## Significance

**The PaperProxy advantage is regularization, not graph-type matching.**

The analytical multinomial formula at low p_eff produces an extremely smooth
N(c';d,c) that gives a well-behaved proxy landscape. The optimizer easily finds
a good optimum in this smooth landscape, and this optimum transfers well to
the true QAOA because:
1. The proxy landscape shape (location of peaks/valleys) is approximately correct
2. The smoothness prevents the optimizer from overfitting to noise

Empirical N(c';d,c), even averaged over 100 instances, retains local structure
that creates minor landscape features the optimizer can get stuck on.

**Implications for the paper:**
- This is actually a POSITIVE result for the proxy approach: the simple
  analytical formula, originally derived for ER, works as a *general-purpose
  regularizer* for any graph family, as long as p_eff is kept low (0.1-0.3).
- The paper should recommend: for non-ER graphs, use PaperProxy with p_eff=0.2
  rather than computing empirical homodist. This is cheaper and gives better results.
- The sampling estimator is still valuable for cases where even PaperProxy fails
  (very dense graphs, p_eff>0.5 regime).

## Next Steps Arising
- [ ] Test PaperProxy with low p_eff (0.1-0.2) at larger n on non-ER. If this
  generalizes, it simplifies the paper's recommendation.
- [ ] Investigate whether the p_eff=0.1-0.2 recommendation also works for
  higher depths (p=5-10).
