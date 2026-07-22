# 024 — The parameter-space findings persist at n = 16

**Question:** E014/E017 (argmax displacement predicts regret), E020
(normalize at depth), and E022 (signed anatomy) were established at
n = 12, 14. Do they persist one size up?

**Answer: Yes, all three, with the same structure (35 instances: 7 families
x 5, n = 16, p = 1 and 3, same seed scheme as E002/E010).**

## Results (results.csv)

1. **Argmax displacement still predicts regret**: pooled Spearman
   rho(regret, displacement) = 0.84 at p=1 and 0.76 at p=3 (normalized
   objective convention). The pooled fidelity-deficit correlation is
   similar here (0.85/0.81), but E017 showed pooled fidelity correlations
   are family-confounded; with 5 instances per cell a within-cell
   decomposition is underpowered at this size, so displacement's
   depth-stable, control-robust status rests on the n = 12-14 analysis.
2. **The normalization rule replicates**: p=3 regret 0.090 (unnormalized)
   vs 0.060 (normalized), normalized better on 32/35; at p=1 the
   unnormalized convention keeps its edge (0.040 vs 0.061, normalized
   better on 2/35). Same directions and similar magnitudes as n = 12-14
   (0.080 -> 0.044).
3. **The signed anatomy replicates**: overprediction at 69/70 measured
   argmax points; winner's-curse ratio |e(theta-hat)|/|e(theta*)| = 4.2
   (p=1) and 9.4 (p=3); margin ~0.04.

## Paper impact

One stability sentence added to the E014/E017 paragraph and one to the
normalization paragraph in section 5.4.

Reproduce: `JULIA_NUM_THREADS=auto julia --project research/experiments/024_n16-stability/run.jl` (~1 h).
