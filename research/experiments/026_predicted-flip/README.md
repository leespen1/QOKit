# 026 — The dispersion theory's prediction, tested

**Question.** E025 says the transfer↔proxy boundary is order-unity weight
dispersion (CV²), not infinite variance. On Pareto, CV² = 1/(α(α−2)) =
2.27 / 1.04 / 0.64 at α = 2.2 / 2.4 / 2.6 — the region E023's grid
skipped. **Pre-registered prediction (written before the run): the flip
sits between α = 2.2 (proxy side) and α = 2.4–2.6 (transfer side).**

**Answer. Confirmed on dense ER — the flip lands at α ≈ 2.4 (CV² ≈ 1.0),
dead center of the pre-registered window; the sparse family's threshold
needs more dispersion (3-regular already flipped before α = 2 per E023),
so the order-unity law carries a family-dependent constant.**

p=1 cell means (array 11756459): ER(0.5) — proxy 0.024 vs sources
0.043–0.068 at α=2.2; a statistical tie 0.028 vs 0.027 at α=2.4; best
source ahead 0.015 vs 0.022 at α=2.6. 3-regular — transfer wins all three
(0.003–0.010 vs proxy ~0.05), consistent with its E023 flip below α=2.
Reading: the dispersion framing holds and is *predictive* (the dense-ER
flip was called before the run); the threshold constant is family-
dependent, larger for sparse regular graphs whose fixed degree gives
concentration more to hold onto.

## Method

E023's harness verbatim, α ∈ {2.2, 2.4, 2.6}, both families, n=16, p=1
and p=3, 10 instances/cell, fresh seeds.

## Reproduce

```bash
cd research/experiments/026_predicted-flip && sbatch run.sb
E26_SMOKE=1 julia --project research/experiments/026_predicted-flip/run.jl  # smoke
```
