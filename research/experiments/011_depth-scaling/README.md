# 011 — E3.2: how does accumulated leakage grow with depth at scale?

**Question:** For p = 30 linear-ramp schedules at n ∈ {16, 18, 20}, does the
accumulated leakage Σλ grow sublinearly in depth at working (small/moderate)
ramps, and how does the growth rate scale with n and graph family?

**Answer: PENDING — running as a 3-task Slurm array on MSU HPCC (one task per
n; CPU only).**

## Method

Extends experiment 003 to larger n and deeper p: per-layer λ_ℓ, Σλ, distance,
overlap, and compressed norm of the exact compressed trajectory
(`compressed_qaoa_trajectory`) for 7 families × 10 instances × 4 ramps
(small/moderate/large/extreme), p = 30. Shared instance seeds.

## Result

*(pending)*

## Reproduce

```
sbatch research/experiments/011_depth-scaling/run.sb        # on HPCC
E1_SMOKE=1 julia --project research/experiments/011_depth-scaling/run.jl    # local smoke
```

Output: `results_task<ID>.csv` per array task (long format, one row per layer).
