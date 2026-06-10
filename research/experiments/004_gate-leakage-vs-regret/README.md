# 004 — GATE: does leakage rank families by parameter-setting regret, and how much regret does the analytical N add?

**Question:** Across graph families, do (i) leakage along the proxy-chosen
schedule, (ii) proxy-state fidelity, and (iii) parameter-setting regret rank
families identically — and how much extra regret does the analytical PaperProxy
N (the only N available without exponential work) add on top of the exact
compression?

**Answer: At p=1 yes — leakage ranks families by regret at Spearman ρ =
0.86–0.96 (gate criterion |ρ| ≳ 0.8 met); at p=3 the family-level signal
vanishes (ρ ≈ 0, but both leakage and regret spreads are narrow along
proxy-chosen schedules, so this is "no discriminating signal," not a
contradiction). Unexpectedly, the analytical PaperProxy beats the exact
compression on every family EXCEPT ER(0.5), where its argmax lands on an
unphysical norm-inflated artifact — which a zero-cost norm filter (a direct
corollary of Theorem 1) repairs.**

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

Family means (n=14 shown; n=12 similar, full table in `results.csv`):

| family        | Σλ (p=3) | regret_emp p1 | regret_paper p1 | regret_emp p3 | regret_paper p3 |
|---------------|----------|---------------|-----------------|---------------|-----------------|
| ER(0.5)       | 0.278    | 0.030         | **0.132**       | 0.080         | 0.117           |
| ER(0.25)      | 0.367    | 0.036         | 0.011           | 0.089         | 0.041           |
| BA(k=2)       | 0.376    | 0.038         | 0.011           | 0.080         | 0.031           |
| BA(k=4)       | 0.260    | 0.028         | 0.013           | 0.087         | 0.066           |
| WS(k=4;b=0.1) | 0.265    | 0.031         | 0.021           | 0.084         | 0.056           |
| WS(k=4;b=0.5) | 0.302    | 0.035         | 0.014           | 0.093         | 0.050           |
| 3-regular     | 0.351    | 0.038         | 0.010           | 0.092         | 0.041           |

1. **Gate verdict.** Family-level Spearman(mean Σλ, mean regret_emp): p=1
   ρ = 0.96 (n=12) and 0.86 (n=14) — **pass**. p=3: ρ ≈ 0 — but along
   proxy-chosen ramps both Σλ (0.26–0.40) and regret (0.06–0.09) barely vary
   across families, so there is no family-level signal to rank rather than an
   anti-correlation. The proxy picks small-angle schedules everywhere, which
   equalizes leakage; the depth/regret question moves to schedule-diverse
   measurements (E2.1) and larger n (E3.1).
2. **H-ER-specific is dead.** PaperProxy with effective edge probability
   p_eff = 2m/(n(n−1)) achieves regret ≈ 0.01–0.02 at p=1 on BA, WS,
   3-regular, and sparse ER — *better* than the exact compression. The smooth
   analytical model acts as a regularizer for parameter choice (consistent
   with the earlier P(c′) investigation). The old log's "proxies fail on
   non-ER graphs" is contradicted on every non-ER family tested.
3. **The ER(0.5) anomaly is an unphysical-artifact story, not a model-quality
   story.** On dense ER(0.5) the analytical-N proxy's grid argmax lands where
   the proxy state's norm has inflated ~7.5× and the predicted ⟨C⟩ = 93
   exceeds the number of edges (38). Theorem 1 makes the exact compression
   contractive — the proxy norm can only decay — so **norm inflation is a
   free certificate of model error**. Filtering grid points to ‖φ‖² ≤ 1.05
   restores PaperProxy on the diagnosed instance to regret ≈ 0.03
   (`diagnose_paper_artifact.jl`). Sparse families have even larger inflation
   regions (norm² up to 617 on 3-regular) but their raw argmaxes happen to
   avoid them — the filter is the principled guard everywhere. Experiment 005
   re-runs the paper arm with this filter across all 420 instances.

## Caveats

- p=3 ceilings/argmaxes use the shared 8⁴ ramp grid (as exp 002).
- Leakage at p=3 is measured along the *proxy-chosen* schedule only; exp 003
  measured fixed schedules and found strong family discrimination — the two
  views answer different questions.
- The PaperProxy arm uses the multinomial formula outside its ER derivation
  domain on BA/WS/3-regular (that is the point of the test).

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/004_gate-leakage-vs-regret/run.jl
```

Seed 20260611 (same instance set as experiment 002). Output: `results.csv`.
Smoke test: prefix `E1_SMOKE=1`.
