# 010 — E3.1: does the Phase-1 picture survive at n = 16–18?

**Question:** At n = 16 and 18 (vs 12–14 in the gate experiment), do the
headroom/regret structure, the near-losslessness of compression at p=1, and
the leakage diagnostics persist — and does the sampled-N proxy (S = 10 per
class, the scalable route) choose parameters as well as the exact-N proxy?

**Answer: PENDING — running as a 14-task Slurm array on MSU HPCC
(one task per (family, n) pair, 20 instances each, A100 ceilings).**

## Method

Mirrors experiment 004 at larger n: per instance, true-QAOA grid ceilings
(p=1 40×40; p=3 ramp 8⁴) on GPU; AR at the argmax of the exact-compression
proxy with exact empirical N (GPU homodist) and with sampled N
(`sampled_homogeneous_distribution`, S=10); uniform/balanced baselines;
Σλ/overlap/distance along the sampled-N-chosen p=3 schedule. Seeds follow the
shared scheme (SEED + 10000·fam + 100·n + inst), so instances extend the
existing family panels.

## Result

*(pending)*

## Reproduce

```
sbatch research/experiments/010_scaleup-ranking/run.sb      # on HPCC
E1_SMOKE=1 julia --project research/experiments/010_scaleup-ranking/run.jl  # local smoke
```

Output: `results_task<ID>.csv` per array task; concatenate for analysis.
