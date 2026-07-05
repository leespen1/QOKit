# 029 — Is the deep high-dispersion corner a grid artifact or intrinsic?

**Question.** E028 found every method bad (regret 0.03–0.3) at (high
weight dispersion × p=10–20 ramps). Coarse-grid artifact, or intrinsic
hardness?

**Answer. Grid artifact.** Compass-search refinement on the true
landscape from ANY grid start — binned proxy, best transfer, universal,
or the grid ceiling itself — converges to essentially one optimum:
continuous ramp endpoints (4 parameters) already cut regret-vs-best-known
from 0.14–0.21 to 0.010–0.037, and freeing all 2p angles finishes the job
(0.000–0.006 from every start, both cells). Even the 8⁴ grid ceiling sat
0.020–0.032 below the refined optimum, so E028's regrets were measured
against a depressed reference. Arrays 11773952 (~34 min/task).

**Reading for the paper:** the deep high-dispersion landscape is not
trap-riddled at these sizes — it is *sharp*: optima fall between the
endpoint-grid points. Wherever true-objective evaluations are affordable
(~2000 statevector calls here), cheap local refinement erases the corner
and the choice of initializer barely matters; the corner is hard only in
the zero-refinement, grid-choice-only setting — e.g. on hardware where
every evaluation costs shots.

## Method

E028's two worst cells (pareto15 × ER(0.5) and logn2.5 × 3-regular, n=16,
p=20, same seeds, 5 instances). Four grid starts per instance; two-stage
compass search (ramp-endpoint 4D, then full-2p 40D; ≤2000 evaluations per
stage); V* = best value found for that instance.

## Reproduce

```bash
cd research/experiments/029_corner-refinement && sbatch run.sb
E29_SMOKE=1 julia --project research/experiments/029_corner-refinement/run.jl  # smoke
```
