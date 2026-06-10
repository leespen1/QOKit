# 010 — E3.1: does the Phase-1 picture survive at n = 16–18?

**Question:** At n = 16 and 18 (vs 12–14 in the gate experiment), do the
headroom/regret structure, the near-losslessness of compression at p=1, and
the leakage diagnostics persist — and does the sampled-N proxy (S = 10 per
class, the scalable route) choose parameters as well as the exact-N proxy?

**Answer: The structure survives scale-up essentially unchanged — regret at
n=16–18 (p=1: 0.03–0.05; p=3 ramps: 0.08–0.11) matches the n=12–14 values
with no growth trend, value-added stays positive at the family level
everywhere, and the sampled-N proxy (S=10) chooses parameters within ~0.01 AR
of the exact-N proxy with systematically slightly *lower* regret (−0.003 to
−0.006 pooled) — the scalable route is validated, and sampling even acts as a
mild regularizer.**

## Method

Mirrors experiment 004 at larger n: per instance, true-QAOA grid ceilings
(p=1 40×40; p=3 ramp 8⁴) on GPU; AR at the argmax of the exact-compression
proxy with exact empirical N (GPU homodist) and with sampled N
(`sampled_homogeneous_distribution`, S=10); uniform/balanced baselines;
Σλ/overlap/distance along the sampled-N-chosen p=3 schedule. Seeds follow the
shared scheme (SEED + 10000·fam + 100·n + inst), so instances extend the
existing family panels.

## Result

280 instances (7 families × n ∈ {16,18} × 20). Family means (regret =
ceiling − AR; exact = instance-exact N via GPU homodist; samp = sampled N,
S=10):

| family        | n  | VA (p=1) | regret p1 exact | regret p1 samp | regret p3 exact | regret p3 samp |
|---------------|----|----------|-----------------|----------------|-----------------|----------------|
| ER(0.5)       | 16 | 0.026    | 0.032           | 0.031          | 0.080           | 0.077          |
| ER(0.5)       | 18 | 0.023    | 0.033           | 0.031          | 0.091           | 0.080          |
| ER(0.25)      | 16 | 0.085    | 0.039           | 0.040          | 0.091           | 0.090          |
| ER(0.25)      | 18 | 0.085    | 0.045           | 0.044          | 0.092           | 0.092          |
| BA(k=2)       | 16 | 0.090    | 0.041           | 0.037          | 0.085           | 0.080          |
| BA(k=2)       | 18 | 0.088    | 0.048           | 0.042          | 0.097           | 0.089          |
| BA(k=4)       | 16 | 0.039    | 0.038           | 0.035          | 0.090           | 0.085          |
| BA(k=4)       | 18 | 0.043    | (rerun pending) | 0.037          | (rerun pending) | 0.078          |
| WS(k=4;b=0.1) | 16 | 0.087    | 0.034           | 0.034          | 0.080           | 0.073          |
| WS(k=4;b=0.1) | 18 | 0.089    | 0.041           | 0.036          | 0.080           | 0.072          |
| WS(k=4;b=0.5) | 16 | 0.084    | 0.044           | 0.042          | 0.094           | 0.084          |
| WS(k=4;b=0.5) | 18 | 0.089    | 0.049           | 0.044          | 0.099           | 0.093          |
| 3-regular     | 16 | 0.127    | 0.046           | 0.042          | 0.105           | 0.092          |
| 3-regular     | 18 | 0.125    | 0.055           | 0.050          | 0.107           | 0.101          |

1. **Regret does not grow with n.** p=1 regret 0.032–0.055 at n=16–18 vs
   0.028–0.038 at n=12–14 (exp 004); p=3 ramps 0.08–0.11 vs 0.06–0.09. The
   compression's parameter-setting quality is scale-stable across the
   accessible range.
2. **Sampled-N is fully validated for parameter setting**: pooled mean
   |AR_samp − AR_exact| = 0.011 (p=1) / 0.016 (p=3) per instance, and pooled
   regret is *lower* for sampled N by 0.003–0.006 — stratified sampling
   smooths instance noise the way the analytical model does, without the
   analytical model's dense-ER artifact. This is the practical recipe the
   filter arc (exps 005–006) failed to find: **use sampled N (S≈10) +
   empirical P.**
3. Value-added stays positive at family level everywhere (0.023–0.127, dense
   families lowest); one ER(0.5) instance posts the campaign's first
   marginally negative instance-level value-added (−0.007) — dense-ER
   headroom is thin and should be stated honestly.
4. **Family regret differences compress at scale** (all families 0.03–0.05
   at p=1), so leakage-vs-regret *rank* correlations are weaker than at
   n=12–14 (ρ ≈ 0.7/0.5) — the same narrow-dynamic-range regime as the
   gate's p=3 cell, not a sign reversal. The ranking claim in the paper
   should be scoped to where regret differences are resolvable.

## Operational notes

- 8 of the original 14 tasks silently fell back to CPU because the @gpu
  environment's CUDA JLL preferences were baked on a driver-less dev node;
  fixed by `CUDA.set_runtime_version!(v"13.0")` in the @gpu env and rerun
  (job 9634813). One rerun task hit a transient concurrent-precompile LLVM
  error; resubmitted solo (job 9636873). The CPU fallback produced valid
  (slower) ceilings — only the exact-N validation columns were affected.

## Reproduce

```
sbatch research/experiments/010_scaleup-ranking/run.sb      # on HPCC
E1_SMOKE=1 julia --project research/experiments/010_scaleup-ranking/run.jl  # local smoke
```

Output: `results_task<ID>.csv` per array task; concatenate for analysis.
