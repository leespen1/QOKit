# 027 — Does the map survive a different problem? Max-3-XOR

**Question.** Every measurement in the paper is MaxCut; Theorems 1–3 are
problem-agnostic and Max-3-XOR is sud2024's second family (integer costs,
native classes). Do the map's pillars replicate: near-lossless p=1
compression, sampled-N ≈ exact-N, and transfer dominance on unweighted
instances?

**Answer.** *Pending — submitted 2026-07-05.*

## Method

Random 3-XOR instances (distinct variables per clause, random parity),
m ∈ {4n, 8n} × n ∈ {14, 16}, 10 instances. (i) One-layer leakage at the
three standard angles. (ii) Parameter setting at p=1 (γ ∈ [0, 2π], 40×40 —
wider γ than MaxCut since the cost operator's spectrum differs) and p=3
ramps (γ ∈ [0.05, 3.2]): exactN, sampledN (S=10), transfer×3 (n=10
same-density sources), universal — against true grid ceilings, regret ≥ 0
by construction.

## Reproduce

```bash
cd research/experiments/027_max3xor && sbatch run.sb
E27_SMOKE=1 julia --project research/experiments/027_max3xor/run.jl  # smoke
```
