# 006 — Does a physicality cap (predicted ⟨C⟩ ≤ m) repair the analytical proxy without collateral damage?

**Question:** Restricting the analytical PaperProxy's argmax to grid points
whose predicted ⟨C⟩ does not exceed the number of edges (an impossibility
bound, unlike experiment 005's calibration-sensitive norm tolerance) — does
that fix the dense-ER(0.5) spurious-peak failure while leaving the
already-good sparse-family argmaxes untouched?

**Answer: No — NEGATIVE RESULT closing the filter arc (004→005→006): on sparse
families even the correctly-located peak carries a predicted ⟨C⟩ above the
physical maximum (the model's values are globally inflated), so the cap rejects
the good peak and regret explodes (0.01 → 0.18–0.25 at p=1), while on dense
ER(0.5) the spurious-peak region extends below m and the cap barely helps
(0.159 → 0.149). The analytical proxy's predicted values carry no usable
absolute information on any family; only its argmax location does.**

## Why this matters

Experiment 005 established that the analytical model's *calibration* (values,
norms) is broken on sparse families even where its *argmax location* is nearly
perfect, so calibration-based filters are catastrophic. A cap at the physical
maximum rejects only predictions that are certainly wrong, never a well-located
peak with merely inflated height — unless the inflation pushes it past m, which
is exactly the ER(0.5) failure mode. If this works, the paper's practical
recommendation becomes: *use the analytical N anywhere, but never trust a
prediction that exceeds the physical maximum.* A tighter heuristic cap
(0.75·m) is included as a sensitivity check; when a cap rejects the whole grid
the selector falls back to the raw argmax (recorded via frac columns).

## Method

Same 420 instances and grids as experiments 004/005. PaperProxy arm only:
AR at argmax under no cap / cap ⟨C⟩ ≤ m / cap ⟨C⟩ ≤ 0.75·m, at p=1 and p=3.

## Result

Mean regret (ceiling − AR), n=12+14 pooled:

| family        | raw p1 | phys p1 | tight p1 | frac_phys p1 | raw p3 | phys p3 | frac_phys p3 |
|---------------|--------|---------|----------|--------------|--------|---------|--------------|
| ER(0.5)       | 0.159  | 0.149   | 0.148    | 0.83         | 0.128  | 0.170   | 0.53         |
| ER(0.25)      | 0.011  | 0.187   | 0.180    | 0.06         | 0.044  | 0.044*  | 0.00         |
| BA(k=2)       | 0.012  | 0.253   | 0.231    | 0.10         | 0.035  | 0.035*  | 0.00         |
| BA(k=4)       | 0.012  | 0.066   | 0.076    | 0.86         | 0.076  | 0.169   | 0.36         |
| WS(k=4;b=0.1) | 0.015  | 0.178   | 0.346    | 0.24         | 0.056  | 0.056*  | 0.00         |
| WS(k=4;b=0.5) | 0.010  | 0.167   | 0.332    | 0.24         | 0.051  | 0.051*  | 0.00         |
| 3-regular     | 0.010  | 0.209   | 0.313    | 0.05         | 0.043  | 0.043*  | 0.00         |

\* frac_phys = 0: every p=3 grid point predicted ⟨C⟩ > m, so the selector fell
back to the raw argmax (no change by construction).

- On sparse families at p=1 only 5–24% of grid points are "physical" — and the
  good argmax is usually NOT among them. The model inflates its predictions
  everywhere, including at the correct peak; rejecting impossible values
  therefore rejects the answer.
- On ER(0.5) most points are physical (83%) including most of the spurious
  region; the cap removes only the artifact's extreme tip.
- Together with experiment 005: **value- and norm-based vetoes cannot repair
  the analytical proxy.** Its usable signal is purely *ordinal/positional*
  (where the peak is), and that signal is excellent off dense graphs and
  corrupted on dense ER(0.5) by a competing spurious peak at large β.

## Where this leaves the practical question

Remaining candidate guards (future experiments, not committed): restricting to
the small-angle trust region where the proxy's derivation holds (the spurious
peak sits at β ≈ 1.17, far outside it — experiment 007's (γ,β) leakage anatomy
will say whether the trust region is definable from leakage itself), or
cross-checking the analytical argmax against a cheap sampled-N proxy
(hypothesis H-sample). For the paper, the filter arc is itself a §6 result:
argmax robustness with a sharply characterized failure mode.

## Caveats

- The tight cap (0.75·m) was a sensitivity check; its added damage on
  WS confirms value-based rejection harms well-located peaks.
- frac_phys at p=3 shows the value decalibration grows with depth (predicted
  values exceed m everywhere on sparse families).

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/006_physicality-filter/run.jl
```

Seed 20260611 (same instance set as 002/004/005). Output: `results.csv`.
Smoke test: prefix `E1_SMOKE=1`.
