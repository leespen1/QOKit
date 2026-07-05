# 026 — The dispersion theory's prediction, tested

**Question.** E025 says the transfer↔proxy boundary is order-unity weight
dispersion (CV²), not infinite variance. On Pareto, CV² = 1/(α(α−2)) =
2.27 / 1.04 / 0.64 at α = 2.2 / 2.4 / 2.6 — the region E023's grid
skipped. **Pre-registered prediction (written before the run): the flip
sits between α = 2.2 (proxy side) and α = 2.4–2.6 (transfer side).**

**Answer.** *Pending — submitted 2026-07-05.*

## Method

E023's harness verbatim, α ∈ {2.2, 2.4, 2.6}, both families, n=16, p=1
and p=3, 10 instances/cell, fresh seeds.

## Reproduce

```bash
cd research/experiments/026_predicted-flip && sbatch run.sb
E26_SMOKE=1 julia --project research/experiments/026_predicted-flip/run.jl  # smoke
```
