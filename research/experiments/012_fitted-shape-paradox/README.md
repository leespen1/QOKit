# 012 — The fitted-shape paradox (E2.4)

**Question:** Does fitting a proxy shape (Triangle/Normal) to the empirical
N(c';d,c) by entrywise MSE actually improve QAOA parameter setting — and if not
(old-log hypothesis H2), does a dynamics-weighted error norm explain why?

**Answer: H2 confirmed, strongly — and our candidate explanation is half wrong.
Entrywise-MSE fitting is harmful (Triangle: mean regret 0.054 → 0.211, worse on
136/140) or useless (Normal: MSE improves 10×, the selected parameters never
change). But the dynamics-weighted *amplitude* error does not rescue the
correlation either: regret is an argmax-transfer quantity, strongly predicted
by the model's argmax displacement (Spearman ≈ 0.7) and by no scalar mismatch
norm we tested (entrywise MSE ≈ −0.1, amplitude error ≈ −0.2, global landscape
correlation ≈ +0.1).**

## Method

For each instance of the seven standard graph families (same seeds as
experiments 002/004; n ∈ {12, 14}, first 10 instances each), at p=1 on the full
40×40 (γ, β) grid:

1. Compute the instance's exact empirical N(c';d,c) (each attained slice sums
   to 2^n) and empirical P(c'). Running the proxy with this N is exact
   compression (experiment 001) — the zero-model-error anchor.
2. Fit a TriangleProxy (4 parameters) and a NormalProxy (3 parameters) to the
   empirical N by minimizing the slice-normalized entrywise MSE — the same
   objective, smart-random-search optimizer, initial values, and bounds as the
   GRIPS Python code (`sendai_opt.fit_proxy_to_real`,
   `run_proxy_study.proxy_init_and_bounds`). Snapshots at 25% and 50% of the
   iteration budget give, together with the unfitted defaults and the final
   fit, a spectrum of models per shape from coarse to converged. The
   analytical PaperProxy (effective edge probability) is a second anchor.
3. For each of the 10 models per instance, record three error measures
   against the *same* instance:
   - **mse_entry** — the fit objective itself (slice-normalized entrywise MSE);
   - **eps_raw / eps_cal** — dynamics-weighted model error: the one-layer
     proxy amplitudes Q₁ under the model N vs under the empirical N,
     ‖Q₁_model − Q₁_emp‖₂, averaged over the grid and at the empirical proxy's
     argmax. `raw` uses the model N as produced (how fitted proxies were
     historically used); `cal` first rescales every N(c';:,:) slice to sum
     2^n — a calibration that requires no instance data;
   - **regret_raw / regret_cal** — grid-ceiling AR minus real-QAOA AR at the
     model's proxy argmax (all models paired with the empirical P, so only the
     N shape varies).

The paradox is confirmed if, across the model spectrum, regret does not
decrease with mse_entry (H2) while it does track the dynamics-weighted error.
The empirical-N rows must reproduce experiment 004's `ar_emp_p1` for the
shared (family, n, instance) cells — a cross-pipeline consistency check.

## Result

140 instances (7 families × n ∈ {12,14} × 10), 10 models each. The empirical-N
rows reproduce experiment 004's `ar_emp_p1` exactly on all 140 shared cells
(cross-pipeline check). Full numbers: `analysis_summary.txt` and
`landscape.csv`; figures `fig_regret_vs_mse.png`,
`fig_regret_vs_eps.png`, `fig_fit_trajectories.png`.

**H2 (does MSE fitting help?) — no, it hurts or does nothing.**

| model | mean MSE (fit objective) | mean regret (raw) |
|---|---|---|
| emp (anchor) | 0 | 0.031 |
| paper | 2.1e-5 | 0.041 |
| tri_init (unfitted) | 2.8e-5 | 0.054 |
| tri_fit | 2.6e-5 | 0.211 |
| norm_init (unfitted) | 1.4e-4 | 0.162 |
| norm_fit | 2.1e-5 | 0.162 |

- Triangle: fitting improves the objective only marginally (median ratio 0.92)
  and raises regret on 136/140 instances (mean +0.157). The unfitted default
  triangle is the *best* shape model.
- Normal: fitting improves the objective 10× (median ratio 0.098) and the
  selected grid point is *identical* to the unfitted one on 140/140 instances
  (Δregret exactly 0; mid-trajectory checkpoints occasionally wobble to other
  argmaxes on 64/140 and return). The Normal-shape landscape's argmax is
  essentially insensitive to its three parameters.
- Both fits converge long before the 1000-iteration budget (the 50-failure
  stop), so the 25%/50% checkpoints nearly always equal the final fit — the
  usable spectrum per shape is effectively {init, fit}.

**Why (which error norm predicts regret?)** Per-instance Spearman ρ of regret
against candidate predictors, across the 9-model spectrum (mean ± sd over 140
instances):

| predictor | ρ (raw) | ρ (calibrated) |
|---|---|---|
| entrywise MSE (fit objective) | −0.09 ± 0.20 | −0.10 ± 0.19 |
| ‖Q₁_model − Q₁_emp‖ at emp argmax | −0.19 ± 0.49 | −0.04 ± 0.19 |
| 1 − Pearson(model landscape, emp landscape) | +0.02 ± 0.26 | +0.10 ± 0.14 |
| argmax displacement ‖(Δγ, Δβ)‖ | **+0.69 ± 0.18** | **+0.74 ± 0.18** |

The amplitude norm is hostage to scale conventions the argmax ignores: the raw
PaperProxy N has slice sums whose per-instance maxima average ~9.2e8 (absolute max 2.3e10;
vs the exact 2^n), giving an
amplitude error of ~1e5 — five orders of magnitude above every other model —
yet the *lowest* non-emp regret (0.041). Conversely the fitted Triangle has
modest amplitude error and the worst regret. Regret is determined by where the
model landscape's argmax lands, and none of the scalar mismatch norms tested
certify argmax transfer. (The displacement correlation is partly mechanical —
a smooth real landscape makes the value gap grow with distance — but that is
the point: the failure mode of fitted shapes is argmax wander, not amplitude
inaccuracy.)

**Calibration (rescaling every slice to 2^n, no instance data needed) is
argmax-neutral at p=1**: it changed the selected grid point on 0/140 instances
for every Triangle/Normal model and 2% for PaperProxy. At p=1 a per-slice
rescale multiplies each |Q₁(c')| by a (γ,β)-independent factor, which evidently
almost never moves the grid argmax. So the historical "raw" usage of fitted
shapes was not a scale bug — and the practical recipe "calibrate the slices" is
harmless but buys nothing at p=1.

**Practical takeaway for §6:** do not fit proxy shapes by entrywise MSE — use
the unfitted default or the analytical N. A useful fitting objective would have
to target the landscape argmax (e.g. match the landscape near its peak), not
the distribution entries.

## Caveats

- p=1 only; the depth behavior of fitted shapes is not tested here.
- All models are paired with the empirical P, which isolates the N shape but
  differs from experiment 004's PaperProxy pairing (binomial P), so the
  `paper` rows here are not directly comparable to 004's `ar_paper_p1`.
- The fit objective normalizes each N(c';:,:) slice independently (the GRIPS
  Python convention), so the fitted shapes are scale-free; the `raw` dynamics
  inherit whatever overall scale the shape's defaults imply.
- The smart-random-search fitter stops after 50 consecutive failures, so the
  25%/50% checkpoints rarely differ from the final fit; conclusions rest on
  the init-vs-fit contrast, not a fine-grained trajectory.
- The argmax-displacement predictor requires the empirical proxy's argmax, so
  it is a diagnosis, not a usable a-priori certificate.

## Reproduce

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/012_fitted-shape-paradox/run.jl
# smoke test (n=10, 2 instances, small grid):
E1_SMOKE=1 julia --project research/experiments/012_fitted-shape-paradox/run.jl
```

Analysis/figures: `julia --project research/experiments/012_fitted-shape-paradox/analyze.jl`
(reads `results.csv`, writes `analysis_summary.txt` and the figures).
Landscape/argmax metrics: `julia --project research/experiments/012_fitted-shape-paradox/landscape_followup.jl`
(rebuilds model Ns from the params column — no refitting — and writes
`landscape.csv` plus the printed displacement correlations in
`landscape_followup.log`).

Paper figure (`fitted_shape_paradox.png`, copied to the paper repo as
`Figures/generated/fitted_shape_paradox.png`):
`julia --project research/experiments/012_fitted-shape-paradox/make_paper_figure.jl`
