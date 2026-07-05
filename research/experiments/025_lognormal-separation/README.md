# 025 — Tail weight or infinite variance? The lognormal separation

**Question.** E023/E024 put the transfer↔proxy crossover at Pareto α ≈ 2,
the infinite-variance point — but Pareto confounds tail weight with
infinite variance. Lognormal(0, σ) weights keep variance finite at every σ
while edge dominance (max-to-sum ratio) grows without bound. Does the
crossover appear anyway (⇒ mechanism is effective single-edge dominance)
or does transfer hold at every σ (⇒ the variance story sharpens)?

**Answer.** *Pending — submitted 2026-07-05.*

## Method

σ ∈ {0.5, 1.0, 1.5, 2.0, 2.5}, both families, n=16, p=1 and p=3,
10 instances/cell — E023's harness verbatim (per-instance mean-weight-
scaled grids; index transfer ≡ rescaled transfer; binned K=64
exact/sampled; 3 weighted sources; pooled universal; regret vs own-grid
ceiling). For calibration: lognormal squared coefficient of variation
e^{σ²}−1 ≈ 0.3, 1.7, 8.5, 54, 518 across the sweep; Pareto α=2 sits at
the divergence of the same quantity.

## Reproduce

```bash
cd research/experiments/025_lognormal-separation && sbatch run.sb
E25_SMOKE=1 julia --project research/experiments/025_lognormal-separation/run.jl  # smoke
```
