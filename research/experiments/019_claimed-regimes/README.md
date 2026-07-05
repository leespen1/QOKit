# 019 — Do the proxy's claimed regimes actually exist?

**Question.** After E013 (no speed advantage at simulable sizes) and E018
(transfer beats the proxy at every tested size), the paper's practical
claims live in the two regimes those experiments did *not* test: high depth
where transfer showed its only strain, and sizes beyond the O(4ⁿ) exact-N
wall where only the analytical/sampled models exist. Does any proxy variant
earn its keep there?

**Answer.** *Pending — Slurm array 11740470 running (submitted 2026-07-04).*

## Method

Two parts, fresh instances (fixed seeds), every method choosing from the
same schedule grid as the true-QAOA grid ceiling (GPU statevector), so
regret ≥ 0 by construction and all methods are directly comparable.

- **Part A (depth, tasks 1–6):** (ER(0.5), ER(0.25), 3-regular) × n ∈
  {14, 16}, 10 instances, p ∈ {10, 20} linear ramps on the standard 8⁴
  endpoint grid. Dense ER at depth is where E018's transfer strained and
  where compression error is largest — if the proxy wins anywhere at
  simulable sizes, it is here.
- **Part B (size, tasks 7–10):** (ER(0.5), 3-regular) × n ∈ {22, 26},
  5 instances, p = 1 (40×40) and p = 3 ramps (8⁴). Exact N is impossible
  (O(4ⁿ)); ceilings remain computable by GPU sweep (~1–8 min each at n=26).

Methods: `ceiling`, `exactN` (Part A only), `sampledN` (S=10 + empirical
P), `paper_binP` (analytical N, native binomial P), `paper_empP`
(analytical N, empirical P), `transfer_1..3` (angles optimized on the true
landscape of an n=10 same-family source; E018's source-seed convention),
`universal` (single family-agnostic schedule: argmax of mean normalized
⟨C⟩ over all 7 families × 3 sources at n=10).

Output: `results_task<ID>.csv`, long format:
`family,n,instance,seed,m,p,method,ar,ceil,regret`.

## Caveats

- Part A depths use ramp-endpoint grids, so "ceiling" is the best *ramp*,
  not the unconstrained p-parameter optimum — consistent with how the
  paper treats depth throughout.
- `universal` is our own pooled-source construction, not a literature
  schedule; it stands in for "no instance knowledge at all."
- n=26 ceilings take minutes per instance on an A100; instance counts are
  5 per cell, so per-cell means carry larger SEs than exps 002/010.

## Reproduce

```bash
cd research/experiments/019_claimed-regimes && sbatch run.sb
E19_SMOKE=1 julia --project research/experiments/019_claimed-regimes/run.jl  # smoke
```
