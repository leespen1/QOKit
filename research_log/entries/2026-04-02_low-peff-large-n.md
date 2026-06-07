# PaperProxy with Low p_eff at Large n and High Depth

**Date:** 2026-04-02
**Paper section:** scalability
**Tags:** paper-proxy, effective-edge-probability, regularization, depth-scaling
**Status:** complete

## Motivation
The regularization finding (p_eff=0.1-0.2 optimal for non-ER at n=12) needs
validation at larger n and higher depths.

## Setup
- **Script:** `julia/paper_figures/low_peff_large_n.jl`
- **Parameters:** n=12,14,16,18; p=1,3,5,8; p_eff=0.2 vs natural
- **Graph types:** BA(m=2), WS(k=4,p=0.3)

## Key Findings

### Low p_eff advantage is inconsistent across n and depth

| Graph | n | p | PP(natural) | PP(0.2) | Winner |
|-------|---|---|-------------|---------|--------|
| WS | 12 | 3 | 0.843 | **0.876** | PP(0.2) by +3.3% |
| WS | 12 | 5 | 0.870 | **0.902** | PP(0.2) by +3.2% |
| WS | 14 | 3 | 0.840 | **0.877** | PP(0.2) by +3.7% |
| WS | 14 | 5 | 0.892 | **0.899** | PP(0.2) by +0.7% |
| BA | 12 | 3 | 0.860 | **0.869** | PP(0.2) by +0.9% |
| BA | 14 | 3 | 0.851 | 0.852 | ~tie |
| WS | 16 | 3 | 0.849 | **0.856** | PP(0.2) by +0.7% |
| BA | 16 | 5 | **0.857** | 0.840 | PP(nat) by +1.6% |
| BA | 14 | 8 | **0.858** | 0.843 | PP(nat) by +1.5% |

PP(0.2) often wins at p=3, but at higher depths (p=5-8) the results are mixed.
At n=18, PP(natural) and PP(0.2) give identical results — the grid search
resolution may not distinguish them at larger n.

### Transfer still dominates both PP variants

Transfer maintains a ~0.03-0.07 advantage over the better PP variant at
all (n, p) combinations. The PP(0.2) improvement doesn't close this gap.

## Significance
- PP(0.2) is a useful heuristic for moderate depths (p=3) on WS graphs
- Not a reliable general recommendation: depends on graph type, n, and p
- Transfer remains the best method for non-ER parameter setting
- Not a strong enough result for a paper recommendation

## Next Steps Arising
None — this line of investigation has diminishing returns. Focus on the
SampN+EmpP headline result on ER instead.
