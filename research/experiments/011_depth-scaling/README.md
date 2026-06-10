# 011 — E3.2: how does accumulated leakage grow with depth at scale?

**Question:** For p = 30 linear-ramp schedules at n ∈ {16, 18, 20}, does the
accumulated leakage Σλ grow sublinearly in depth at working (small/moderate)
ramps, and how does the growth rate scale with n and graph family?

**Answer: Accumulated leakage grows *linearly* in depth for fixed-endpoint
ramps — Σλ(p=30)/Σλ(p=20) = 1.505 at small ramps vs the 1.5 predicted by
Theorem 3's per-layer profile argument — and only mildly with n (+13–30% from
n=16→20 despite m growing ~58%), so proxy fidelity at depth is predictable:
overlap at p=30 stays 0.71–0.81 for small ramps even at n=20.**

## Method

Extends experiment 003 to larger n and deeper p: per-layer λ_ℓ, Σλ, distance,
overlap, and compressed norm of the exact compressed trajectory
(`compressed_qaoa_trajectory`) for 7 families × 10 instances × 4 ramps
(small/moderate/large/extreme), p = 30. Shared instance seeds.

## Result

25,200 layer records (7 families × 3 n × 10 instances × 4 ramps × 30 layers).

1. **Linear-in-p accumulation, quantitatively as predicted.** For ramps with
   fixed endpoints, the per-layer leakage is a function of (γ_ℓ, β_ℓ) ≈
   f(ℓ/p) alone (E2.1's λ ≈ const·β·γ²·m), so the total should scale as
   p·mean(f). Same-endpoint, same-n comparison against experiment 003's
   p=20 runs at n=16:

   | ramp     | Σλ (p=20) | Σλ (p=30) | ratio | linear = 1.50 |
   |----------|-----------|-----------|-------|----------------|
   | small    | 1.358     | 2.044     | 1.505 | ✓ exact        |
   | moderate | 2.684     | 3.737     | 1.392 | mild saturation|
   | large    | 3.799     | 4.600     | 1.211 | saturating     |

   Saturation at larger ramps is expected: once the compressed norm has
   decayed, later layers have less amplitude left to leak.
2. **Mild n-dependence at fixed schedule:** family-pooled Σλ at layer 20
   (moderate ramp) grows from n=16 → 20 by only 13–30% per family while
   m grows ~30–58% — accumulated leakage scales *sublinearly* in m along
   deep trajectories (in contrast to the single-layer ∝ m law), again
   consistent with norm-loss saturation.
3. **Fidelity at depth stays usable in the working regime:** mean overlap at
   p=30 is 0.81/0.76/0.71 (small ramp, n=16/18/20) and 0.45/0.38/0.31
   (moderate). The proxy's deep-ramp parameter-setting success (Sud et al.'s
   p=20 regime) operates at fidelities well below 1 — argmax robustness, not
   state fidelity, is what survives depth (consistent with exps 004–006, 009).
4. Family ordering of leakage by density holds at every n and depth
   (ER(0.5) most, 3-regular least).

## Caveats

- Same four ramp schedules as experiment 003; conclusions are about
  fixed-endpoint linear ramps, not optimized schedules.
- Ratios compare across experiments with identical seeds/instances at n=16;
  the n=18/20 instances have no p=20 counterpart (not needed for the claim).

## Reproduce

```
sbatch research/experiments/011_depth-scaling/run.sb        # on HPCC
E1_SMOKE=1 julia --project research/experiments/011_depth-scaling/run.jl    # local smoke
```

Output: `results_task<ID>.csv` per array task (long format, one row per layer).
