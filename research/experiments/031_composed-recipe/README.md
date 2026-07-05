# 031 — The composed recipe, end to end

**Question.** The paper's closing advice composes three findings: init
from transfer/universal where family concentration holds (E018/E019) or
from the binned proxy under high dispersion (E022–E026); polish on the
surrogate within the ramp family (E030); confirm with one true
evaluation. Does the fully deployable composition — zero target
evaluations except the confirmation — actually deliver in each regime?

**Answer. Yes, with a two-evaluation guard.** Unweighted (transfer's
home): inits already sit at the ceiling (0.897) and polish *degrades*
them (0.837) — the guard (evaluate init and polished, keep the better)
recovers fully. Proxy's niche (Pareto p=3): polish +0.05 (0.76–0.77 →
0.82). Hard corner (lognormal p=20): polish +0.09–0.12 (0.77–0.80 →
0.90), best-of-two 0.92–0.95. The deployable recipe — map-init,
surrogate ramp-polish, two true evaluations, keep the better — never
hurts and gains exactly where the map says the surrogate has signal.
Array 11782218, all COMPLETED (~4–6 min/task).

## Method

Three regimes: unweighted ER(0.5) p=3 (transfer's home), Pareto(1.5)
ER(0.5) p=3 (the proxy's niche), lognormal(2.5) 3-regular p=20 (the hard
corner); n=16 throughout. Inits: `transfer_1` (single weighted source, no
target evaluations) and `universal` (pooled). Each init polished by
compass search of the 4 ramp endpoints on the binned-sampled proxy's
landscape (K=64, S=10); one true evaluation per stage; reported against
the 8⁴ grid ceiling.

## Reproduce

```bash
cd research/experiments/031_composed-recipe && sbatch run.sb
E31_SMOKE=1 julia --project research/experiments/031_composed-recipe/run.jl  # smoke
```
