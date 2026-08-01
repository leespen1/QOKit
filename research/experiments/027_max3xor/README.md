# 027 — Does the map survive a different problem? Max-3-XOR

**Question.** Every measurement in the paper is MaxCut; Theorems 1–3 are
problem-agnostic and Max-3-XOR is sud2024's second family (integer costs,
native classes). Do the map's pillars replicate: near-lossless p=1
compression, sampled-N ≈ exact-N, and transfer dominance on unweighted
instances?

**Answer. Everything replicates.** Exact-compression p=1 regret
0.031–0.040 (MaxCut magnitude); sampled-N matches (within 0.006) or beats
exact in 7 of 8 cells, strictly beats in 5 (clearest: 8n, n=16, p=3:
0.044 vs 0.086 — the regularization effect's
fourth independent occurrence); transfer and even the pooled universal
schedule are near-perfect on unweighted 3-XOR (regret ≈ 0.000 in most
cells — concentration stronger than MaxCut's); small-angle leakage
doubles when clause density doubles (0.045→0.092 at n=14: the λ ∝ m
density law again). Large-angle leakage is much higher than MaxCut's
(0.84–0.95 at γ=1.0), consistent with 3-local phases lacking the MaxCut
cancellation — the one genuinely MaxCut-specific piece, as the theory
says. Array 11757023, all COMPLETED.

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
