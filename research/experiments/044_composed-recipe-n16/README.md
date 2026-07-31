# E044 — Does the composed recipe (sampled N + normalized) survive n=16?

**Question**: Does E041's headline (sampled S=10 + normalized objective is
the best proxy combination at both depths, beating even exact+normalized)
replicate one size up, at n=16, including the S=3 cell E042 found strongest?

**Answer**: **Yes, emphatically: at n=16 sampled+normalized beats
exact+normalized on 34/35 instances at BOTH depths (sign p=1e-9), the
advantage grows with size (p=1 mean regret 0.018 vs 0.061; p=3 0.024 vs
0.060), and S=3 is stronger still (0.0067 / 0.0211).**

## Method

7 families × n=16 × 5 instances (E041 seed scheme). Variants
{exact, samp10, samp3} × {raw, normalized}; p=1 grid (40×40) and p=3
linear-ramp grid (8⁴); regret vs the true statevector grid ceiling.
Mirrors E040's stability role for the earlier findings.

## Result

Mean regret over 35 instances:

| variant      | p=1 raw | p=1 norm | p=3 raw | p=3 norm |
|--------------|---------|----------|---------|----------|
| exact        | 0.040   | 0.061    | 0.090   | 0.060    |
| sampled S=10 | 0.045   | **0.018**| 0.091   | **0.024**|
| sampled S=3  | 0.046   | 0.0067   | 0.097   | 0.021    |

The exact-N normalized cell reproduces E040 (0.060). Every qualitative
feature of E041/E042 replicates one size up: normalization hurts the exact
compression at p=1 but flips to a large win under sampling; raw-objective
regret is insensitive to sampling; smaller S helps more. The n-trend is
notable: the composed recipe's advantage over exact+normalized widens from
n=12/14 (0.021 vs 0.047) to n=16 (0.018 vs 0.061).

## Caveats

- 5 instances per family (35 total), same as E040's n=16 budget.
- Single sampling draw per (instance, S); E042 quantified draw variance.

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/044_composed-recipe-n16/run.jl   # ~40 min
E44_SMOKE=1 julia --project research/experiments/044_composed-recipe-n16/run.jl              # quick check
```
