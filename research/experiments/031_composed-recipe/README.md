# 031 — The composed recipe, end to end

**Question.** The paper's closing advice composes three findings: init
from transfer/universal where family concentration holds (E018/E019) or
from the binned proxy under high dispersion (E022–E026); polish on the
surrogate within the ramp family (E030); confirm with one true
evaluation. Does the fully deployable composition — zero target
evaluations except the confirmation — actually deliver in each regime?

**Answer.** *Pending — submitted 2026-07-05.*

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
