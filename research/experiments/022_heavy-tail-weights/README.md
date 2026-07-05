# 022 — Strong heterogeneity, fair transfer: the practitioner finish line

**Question.** E021 showed i.i.d. U[0,1] weights erode but don't break
transfer. With genuinely heavy-tailed weights (Pareto α=1.5, infinite
variance; Exp(1) as the milder rung) — and transfer given its
practitioner's fix for free (grid-index transfer over per-instance
mean-weight-scaled γ-grids ≡ rescaled angle transfer) — does the binned
proxy finally win a regime, or does the practitioner section close with
"always transfer"?

**Answer.** *Pending — submitted 2026-07-05.*

## Method

(ER(0.5), 3-regular) × n ∈ {14, 16} × law ∈ {Exp(1), Pareto(1.5)},
10 instances each. Every instance's schedule grid uses γ = γ̂/w̄ (shared
normalized γ̂-grid, its own mean weight w̄), so all methods — ceiling,
binned exact/sampled proxy (K=64, S=10), transfer×3 (weighted n=10
same-family sources on their own scaled grids), pooled universal —
choose from the same per-instance grid, and index transfer implements
mean-weight rescaling exactly. Regret ≥ 0 by construction (asserted).

## Caveats

- Mean-weight rescaling is one practitioner heuristic; total-weight or
  quantile-based rescalings could differ for Pareto (the mean is
  tail-sensitive at α=1.5).
- Pareto instances have occasional dominant edges; the quantile bins
  absorb them by construction, transfer sees them only through w̄.

## Reproduce

```bash
cd research/experiments/022_heavy-tail-weights && sbatch run.sb
E22_SMOKE=1 julia --project research/experiments/022_heavy-tail-weights/run.jl  # smoke
```
