# 004 — GATE: does leakage rank families by parameter-setting regret, and how much regret does the analytical N add?

**Question:** Across graph families, do (i) leakage along the proxy-chosen
schedule, (ii) proxy-state fidelity, and (iii) parameter-setting regret rank
families identically — and how much extra regret does the analytical PaperProxy
N (the only N available without exponential work) add on top of the exact
compression?

**Answer: PENDING — full run queued.**

## Why this matters

This is the Phase-1 **go/no-go gate** (program: gate ~Jul 8). The paper's
central claim is that compression error (leakage) explains when the homogeneous
heuristic works. Exp 003 already gives leakage ⇔ *fidelity* at ρ ≈ 1; the
missing and riskiest link is leakage ⇔ *regret* (a bad-fidelity proxy could
still pick good parameters). Pre-committed criterion: family-level
|Spearman ρ| ≳ 0.8 between mean Σλ and mean regret_emp → proceed theory-led;
otherwise pivot to the fidelity-only framing.

The PaperProxy arm (analytical Binomial/Multinomial N with effective edge
probability p_eff = 2m/(n(n−1)), binomial P) tests exp 002's prediction that
non-ER failure is **model error**: regret_paper − regret_emp is the model
error's cost in AR, measured per family.

## Method

Same 420 instances as experiment 002 (shared seeds). Per instance, at p=1
(40×40 grid) and p=3 (8⁴ ramp grid): grid ceiling, real AR at the empirical-N
proxy argmax (→ regret_emp), real AR at the PaperProxy-N argmax
(→ regret_paper), and Σλ / overlap / distance along the empirical proxy's
chosen schedule.

## Result

*(pending)*

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/004_gate-leakage-vs-regret/run.jl
```

Seed 20260611 (same instance set as experiment 002). Output: `results.csv`.
Smoke test: prefix `E1_SMOKE=1`.
