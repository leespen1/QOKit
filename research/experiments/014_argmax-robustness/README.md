# 014 — At depth, is regret a fidelity problem or an argmax-transfer problem?

**Question:** Exp 004 found that at p=3 the proxy-chosen schedules equalize
leakage across families (leakage⇔regret ρ≈0), so fidelity can no longer be what
sets regret there. Is depth regret governed instead by *argmax transfer* (how far
the proxy's chosen parameters sit from the true optimum) or by *landscape
robustness* (how flat the true peak is)? And does this explain why a worse-fidelity
model (the analytical N) can pick better parameters than the exact compression?

**Answer: Parameter-setting regret is an argmax-transfer (parameter-space)
quantity, decoupled from state fidelity at depth. Pooled over 140 instances,
regret correlates with argmax displacement at both depths (Spearman ρ = 0.76 at
p=1, 0.63 at p=3) but with the fidelity deficit only at p=1 (ρ = 0.39 → −0.02 at
p=3 — fidelity stops predicting regret entirely). The "landscape robustness
governs depth regret" hypothesis is refuted (ρ ≈ −0.18/+0.10, weak and
sign-inconsistent). So the leakage calculus bounds the STATE error, but
parameter-setting quality is a separate, parameter-space matter — exactly why a
low-fidelity model with a well-placed argmax (the analytical N off-ER) beats the
high-fidelity exact compression.**

## Why this matters

This resolves the open framing question in STATUS (how §5 should narrate depth)
and gives §5.4's claim — "parameter setting is an argmax-transfer problem, not a
state-approximation problem" — direct instance-level evidence at depth. It closes
the loop on the 004→005→006→012 arc: across model classes (012) *and* across
depth (014), the only thing that tracks regret is argmax displacement; fidelity
and every value-/norm-based quantity fail.

## Method

The exp 002/004 instance set (same seeds), subset 7 families × n∈{12,14} × 10
instances. For each instance we recompute the true-landscape geometry exp 004 did
not save, at p=1 (40×40 grid over γ∈[0,π], β∈[0,π/2]) and p=3 (8⁴ linear-ramp
endpoint grid):

- **ceiling, ar_emp, regret** — best real-QAOA AR on the grid; real AR at the
  exact-compression proxy's argmax; their difference. (Reproduces exp 004's
  regret_emp as a cross-check.)
- **argmax displacement** — normalized angle distance between the proxy's argmax
  and the true argmax (each angle scaled to its grid range; Euclidean / √(#angles)).
- **robustness (flat-peak fraction)** — fraction of grid points within 0.01 AR of
  the ceiling: large = flat, forgiving peak.
- **overlap** — proxy-state fidelity |⟨ψ|φ⟩|² at the chosen angles (the fidelity
  axis), from `compressed_qaoa_trajectory`.

`analyze.jl` then correlates regret against each, at p=1 and p=3 separately.

## Result

Pooled Spearman ρ(regret, ·) over 140 instances per depth:

| predictor | p=1 | p=3 |
|---|---|---|
| fidelity deficit (1−overlap) | **0.39** | **−0.02**  ← decoupled |
| argmax displacement | **0.76** | **0.63** |
| landscape robustness (flat-peak %) | −0.18 | +0.10  (refuted) |

`argmax_vs_fidelity.png`: (a) regret vs fidelity deficit — a slope at p=1 that
flattens to a high cloud at p=3; (b) regret vs argmax displacement — a positive
trend at both depths.

The reading: as depth grows, the proxy-chosen schedules push fidelity into a
narrow band (CV(overlap) = 0.014 at p=3) that no longer discriminates good from
bad parameter choices, while the *parameter-space* gap between the proxy's argmax
and the truth keeps setting regret. Parameter setting lives in parameter space.

## Caveats

- Subset (10 instances/family, n=12,14); p=1 full 40×40 grid, p=3 the same 8⁴
  ramp grid as exp 002/004.
- The flat-peak robustness measure is coarse at p=3 (≈20 of 4096 grid points fall
  within ε=0.01 of the ceiling), so its null result is weakly resolved; the firm
  conclusions are the fidelity decoupling and the displacement persistence.
- regret↔displacement is partly mechanical (a displaced argmax on a curved peak
  costs regret); the *finding* is the contrast — fidelity decouples while
  displacement persists — not the displacement correlation alone.
- Family-mean fidelity⇔regret at p=3 (ρ=0.61) is confounded by pooling n=12 and
  n=14 (both regret and infidelity rise with n); the clean, unconfounded signal is
  the instance-level pooled ρ=−0.02.

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/014_argmax-robustness/run.jl
julia --project research/experiments/014_argmax-robustness/analyze.jl
julia --project research/experiments/014_argmax-robustness/make_figure.jl
```

Instances seeded by `20260611 + 10000·fam_idx + 100·n + inst` (identical to exp
002/004). Smoke test: `E14_SMOKE=1 julia --project .../run.jl`.
