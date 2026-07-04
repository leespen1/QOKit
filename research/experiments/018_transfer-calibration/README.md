# 018 — Transfer calibration: is the proxy better than just borrowing angles?

**Question.** How does the proxy's parameter-setting regret compare to the
standard cheap alternative — transferring angles grid-optimized on a small
(n=10) instance of the same family? (Paper-critique item T2.4: the one
calibration number a referee will predictably ask for.)

**Answer.** *Pending — Slurm job 11677649 running (submitted 2026-07-04).*

## Method

For every target instance already measured in exps 002 (n=12,14; 30
instances/cell) and 010 (n=16,18; 20 instances/cell), across the seven
standard families: regenerate the instance exactly from the seed recorded in
those experiments' CSVs (edge count asserted against the CSV), then evaluate
real QAOA at angles grid-optimized on the *true* landscape of an n=10 source
instance of the same family — 3 independent sources per family — at p=1
(40×40 grid) and p=3 linear ramps (8⁴ endpoint grid), the same grids as
002/010. Transfer regret is measured against the target's own recorded grid
ceiling, so it is directly comparable, cell by cell, to the exact-compression
and sampled-N proxy regrets already in those CSVs. Because the source angles
live on the same schedule grids as the ceilings, transfer regret is
nonnegative by construction (asserted).

Output: `results.csv`, long format — one row per (target, source, p):
`family,n,instance,seed,source_rep,p,ar_transfer,ar_ceiling,regret_transfer`.

## Caveats

- Sources are n=10 with the *true* (statevector) landscape — the most
  favorable version of transfer (no surrogate error on the source side); this
  makes the comparison conservative *against* the proxy.
- Transfer is within-family with unscaled angles; rescaling rules
  (e.g. γ ∝ 1/√m) could improve the dense→dense transfers and are cited, not
  tested, here.
- p=3 transfers share the ramp-endpoint grid with the ceilings; finer
  schedules could shift both sides.

## Reproduce

```bash
cd research/experiments/018_transfer-calibration && sbatch run.sb
E18_SMOKE=1 julia --project research/experiments/018_transfer-calibration/run.jl  # smoke
```
