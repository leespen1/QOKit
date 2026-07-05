# 022 — Strong heterogeneity, fair transfer: the practitioner finish line

**Question.** E021 showed i.i.d. U[0,1] weights erode but don't break
transfer. With genuinely heavy-tailed weights (Pareto α=1.5, infinite
variance; Exp(1) as the milder rung) — and transfer given its
practitioner's fix for free (grid-index transfer over per-instance
mean-weight-scaled γ-grids ≡ rescaled angle transfer) — does the binned
proxy finally win a regime, or does the practitioner section close with
"always transfer"?

**Answer. The crossover exists: under Pareto(α=1.5) weights the binned
proxy wins every p=1 cell by 2–7× (e.g. 0.012 vs 0.087–0.111 on G(16,0.5))
and is the consistent choice at p=3, where single transfer sources become a
lottery (0.02–0.14). Exponential weights remain transfer's territory.**

Array 11751040, all 8 tasks COMPLETED (~1–3 min each). Cell means in
analysis below. Reading:

- **Exp(1):** transfer still wins everywhere (regret 0.002–0.015 vs binned
  0.028–0.117) — mild i.i.d. heterogeneity behaves like E021's U[0,1].
- **Pareto(1.5), p=1:** binned proxy 0.012–0.022; every transfer source
  0.031–0.099; universal 0.032–0.111. Mean-weight rescaling fails because
  the mean is tail-dominated — one or two giant edges reshape the
  landscape in a way no cross-instance schedule can track, while the
  quantile bins absorb them per instance.
- **Pareto(1.5), p=3:** best-of-3 transfer still edges the binned proxy in
  3 of 4 cells, but with enormous source variance (same-family sources
  range 0.021 → 0.140); the binned proxy sits stably at 0.046–0.067.
  Risk-adjusted, the per-instance method is the defensible choice.
- binned_sampled ≥ binned_exact in most Pareto cells (sampling
  regularizes; now a third independent occurrence).

**The practitioner map is complete:** family concentration ⇒ transfer;
heavy-tailed/heterogeneous weights ⇒ per-instance binned proxy; the
calculus tells you which regime you are in (leakage and bin spread are
measurable before committing).

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
