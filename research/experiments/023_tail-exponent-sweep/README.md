# 023 — Where is the crossover? A tail-exponent sweep

**Question.** E022 found the transfer→proxy ranking inversion at
Pareto(α=1.5) and its absence at Exp(1). Is α=1.5 special, or does the
inversion set in smoothly as tails get heavier — and where?

**Answer. The crossover sits at the infinite-variance boundary α ≈ 2.**
At p=1 (binned-sampled proxy vs best-of-3 rescaled transfer): α=1.2 —
proxy 0.017–0.020 vs 0.035/0.128 (single sources as bad as 0.21); α=1.5 —
proxy 0.012–0.014 vs 0.020/0.086; α=2.0 — 3-regular flips to transfer
(0.007 vs 0.034) while dense ER is a statistical tie (0.028 vs 0.028);
α=3, 5 — transfer wins cleanly (0.0006–0.006 vs 0.029–0.046). The proxy's
regret is nearly α-independent (~0.02–0.03 at p=1) — it is transfer that
moves through it as the tail lightens. At p=3 the same pattern shifts one
rung heavier (proxy wins only at α=1.2; the source lottery is extreme:
same-family sources at α=1.2 range 0.035–0.212). Array 11751921, all 10
tasks COMPLETED.

Reading: for α > 2 the weight distribution has finite variance, instance
landscapes concentrate around the family mean, and borrowing angles works;
at α ≤ 2 single edges dominate instances individually and only a
per-instance method tracks them. The boundary is principled, not tuned.

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
