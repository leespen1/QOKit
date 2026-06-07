# Random Partition Baseline Approximation Ratios

**Date:** 2026-04-02
**Script:** `julia/paper_figures/random_baseline.jl`
**Status:** Complete

## Question

What is the random partition baseline for all experimental configurations?
This contextualizes all prior QAOA improvements — if the random baseline is
already high, QAOA "improvements" operate in a narrow band.

## Setup

- Graph types: ER(0.5), BA(m=2), WS(k=4, p=0.3), ER(p_edge={0.1,0.2,0.3,0.7})
- n = 12-22, 20 instances per configuration
- Unbiased baseline: mean(costs)/max(costs)
- Biased baseline: best over q in {0.1, 0.15, ..., 0.9} with 10K samples each

## Key Results

### Standard graph types (n=12-22)

| Graph | n=12 | n=14 | n=16 | n=18 | n=20 | n=22 |
|-------|------|------|------|------|------|------|
| ER(0.5) | 0.692 | 0.695 | 0.706 | 0.729 | 0.736 | 0.753 |
| BA(m=2) | 0.634 | 0.631 | 0.613 | 0.615 | 0.624 | — |
| WS(k=4) | 0.639 | 0.643 | 0.627 | 0.631 | 0.627 | — |

### ER with varying p_edge (n=18)

| p_edge | Random baseline | Prior QAOA context |
|--------|----------------|-------------------|
| 0.1 | 0.529 | Very low — large headroom for QAOA |
| 0.2 | 0.602 | Low — good target for proxy |
| 0.3 | 0.647 | Moderate |
| 0.5 | 0.731 | High — small improvement band |
| 0.7 | 0.795 | Very high — tiny improvement band |

### Contextualizing Prior Headline Results

At n=22, p=3 on ER(0.5):
- Random baseline: **0.753**
- Transfer: 0.844 → gap above random = **0.091**
- SampN+EmpP: 0.881 → gap above random = **0.128**
- SampN+EmpP's improvement over Transfer (0.037) represents **41% more
  improvement** above random than Transfer alone

At n=18, p=3 on ER(0.5):
- Random baseline: **0.729**
- Transfer: 0.846 → gap = 0.117
- SampN+EmpP: 0.865 → gap = 0.136

## Key Insights

1. **Biased partitions don't help**: For all configurations tested, the best
   biased partition matches the unbiased baseline. The random baseline is
   mean(costs)/max(costs) ≈ m/2/c_opt regardless of bias.

2. **ER(0.5) baselines are high**: At n=22, random already gives 0.753. All
   QAOA methods operate in the 0.75-0.90 range. This means percentage-point
   improvements are actually substantial relative to the available headroom.

3. **Sparse/non-ER baselines are much lower**: BA(m=2) ~0.62, WS(k=4) ~0.63,
   ER(0.1) ~0.52. These have much more headroom for QAOA improvement. A proxy
   method achieving even 0.75 on sparse BA would be very significant.

4. **ER(0.1-0.2) is the most promising testbed**: Random baseline only 0.51-0.60,
   PaperProxy formula supports arbitrary p_edge, and the improvement headroom
   is large. These should be prioritized for the next experiment.

5. **Random baseline should appear on ALL paper figures**: It provides essential
   context. Without it, a reader cannot assess whether 0.88 vs 0.84 is meaningful.

## Implications for Paper

- All approximation ratio figures must include a horizontal line showing the
  random baseline for the configuration being plotted.
- The "real improvement" metric should be reported: (method - random) / (1 - random).
- Sparse ER and non-ER graphs are more interesting testbeds because the random
  baseline is lower, making any proxy improvement more impactful.
