# 024 — Does the heavy-tail crossover survive scale?

**Question.** E022/E023 located the transfer↔proxy inversion at the
infinite-variance boundary (α ≈ 2) at n=14–16. Does the α<2 proxy win
persist at n=22 and 26 — sizes where the exact binned homodist is
O(4ⁿ)-impossible and the proxy side is necessarily the *sampled* binned
variant (K=64, S=10), i.e. the practical method whose niche is claimed?

**Answer.** *Pending — submitted 2026-07-05.*

## Method

(ER(0.5), 3-regular) × n ∈ {22, 26} × Pareto α ∈ {1.5, 3.0}, 5 instances,
p=1 on per-instance mean-weight-scaled grids (γ̂ ∈ [0, 2π], 40×40).
Ceilings by GPU statevector sweep; methods: `binned_sampled` (the only
feasible proxy at this size), `transfer_1..3` (weighted n=10 sources,
index transfer ≡ rescaled transfer), `universal`. Regret ≥ 0 by
construction.

## Caveats

- 5 instances/cell (GPU-time bound): read cell means with ±SE in mind.
- p=1 only; the depth dimension at scale remains untested for heavy tails.

## Reproduce

```bash
cd research/experiments/024_heavy-tail-at-scale && sbatch run.sb
E24_SMOKE=1 julia --project research/experiments/024_heavy-tail-at-scale/run.jl  # smoke
```
