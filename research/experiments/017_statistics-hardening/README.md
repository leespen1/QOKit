# 017 — Statistics hardening: do the headline correlations survive uncertainty quantification and family controls?

**Question:** The journal-readiness review (2026-07-22) flagged that the
paper's ~10 reported Spearman coefficients carry no confidence intervals, and
that the pooled 140-instance correlations (E014) and 1260-row correlations
(E012) conflate between-family and within-family variation (a Simpson's-type
confound). Do the paper's claims survive bootstrap CIs, cluster bootstrap
over family-size cells, within-cell correlations, and cell demeaning?

**Answer: The argmax-displacement claim survives every control; the fidelity
story needs (and now has) a more careful statement; everything else holds
with quantified uncertainty.**

## Results (results.csv; B = 5000 bootstrap resamples, seed 20260722)

1. **Argmax displacement predicts regret, robustly, at both depths.**
   Pooled rho = 0.74 (95% CI 0.64-0.82) at p=1 and 0.65 (0.54-0.74) at p=3;
   cluster bootstrap over the 14 family-size cells leaves the CIs off zero
   (0.56-0.87, 0.52-0.74); mean within-cell rho 0.65 and 0.56; cell-demeaned
   pooled rho 0.63 and 0.67. This is the paper's most robust empirical claim.
2. **The fidelity correlations are family-confounded, in both directions.**
   At p=1 the pooled rho = 0.39 (0.24-0.53) collapses to 0.07 (n.s.) within
   cells: fidelity's apparent p=1 predictiveness is mostly family identity.
   At p=3 the pooled rho = -0.02 hides a moderate within-cell correlation of
   0.46 (0.21-0.72): proxy-chosen p=3 schedules compress fidelity into a
   narrow band ACROSS families (as the paper's mechanism says), which zeroes
   the pooled statistic, but within a family-size cell fidelity still carries
   some signal. The honest summary, now in the paper: fidelity is not a
   stable predictor at either depth once family identity is controlled;
   argmax displacement is, at both.
3. **Landscape flat-peak robustness predicts nothing** at either depth under
   any control (|rho| <= 0.24, all CIs straddle or nearly straddle 0).
4. **E012 norms:** entrywise MSE rho vs regret is ~0.03 pooled over 1260
   (instance, model) rows, ~0.02 within cells, 0.004 demeaned: no predictive
   power, robustly. Amplitude error at the argmax is consistently mildly
   NEGATIVE (-0.26 pooled, -0.21 within, -0.30 demeaned): bigger amplitude
   error goes with LOWER regret across the model spectrum, because the
   best-argmax model (analytical) has the worst calibration. Reinforces the
   paper's "no norm certifies argmax transfer."
5. **E004 7-point ranking now carries CIs:** rho = 0.96 with instance-level
   bootstrap 95% CI [0.64, 0.96] at n=12; 0.86 with [0.61, 0.93] at n=14.
   Both exclude 0 decisively; the paper additionally notes the rho >= 0.8
   pre-commitment.
6. **Depth-accumulation ratio:** per-instance Sum-lambda(p=30)/Sum-lambda(p=20)
   over the 70 matched n=16 small-ramp runs (E011 vs E003, identical seeds)
   is 1.536 +- 0.262 (SD). The paper's 1.505 is the ratio of pooled sums;
   both are consistent with the predicted 1.5.
7. **E010 sign test:** sampled-N regret is lower than exact-N in 26 of 28
   family-size cells; one-sided binomial p = 1.5e-6 treating cells as
   independent (they share ceiling-computation regimes, so this is
   indicative, not exact).

## Method

Pure-Julia analysis over the committed CSVs of E003/E004/E010/E011/E012/E014
(no new simulation). Tie-averaged Spearman; percentile bootstrap
(instance-level and cluster-level, resampling whole family-size cells);
within-cell rho = mean over the 14 cells of 10 instances (cells with zero
variance in a variable are skipped: 1 such cell for argmax_disp at p=1);
cell-demeaned rho = Spearman after subtracting cell means from both
variables.

Reproduce: `julia --project research/experiments/017_statistics-hardening/run.jl`

## Paper impact

- Section 5.4 E014 paragraph rewritten with CIs and the within-family
  decomposition (both confound directions stated).
- Section 5.3 ranking sentence now carries the bootstrap CIs.
- The 1.505 sentence now carries the per-instance spread.
- The 26/28 claim now carries the sign-test p-value with the independence
  caveat.
