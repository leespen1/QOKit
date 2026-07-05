# 023 — Where is the crossover? A tail-exponent sweep

**Question.** E022 found the transfer→proxy ranking inversion at
Pareto(α=1.5) and its absence at Exp(1). Is α=1.5 special, or does the
inversion set in smoothly as tails get heavier — and where?

**Answer.** *Pending — submitted 2026-07-05.*

## Method

Pareto α ∈ {1.2, 1.5, 2.0, 3.0, 5.0} (mean α/(α−1): 6, 3, 2, 1.5, 1.25;
variance infinite for α ≤ 2), both families at n=16, 10 instances, p=1 and
p=3 — otherwise identical to E022 (per-instance mean-weight-scaled grids;
index transfer = the practitioner's rescaling; binned K=64 exact/sampled;
3 weighted sources; pooled universal; regret vs each instance's own grid
ceiling, nonnegative by construction).

## Reproduce

```bash
cd research/experiments/023_tail-exponent-sweep && sbatch run.sb
E23_SMOKE=1 julia --project research/experiments/023_tail-exponent-sweep/run.jl  # smoke
```
