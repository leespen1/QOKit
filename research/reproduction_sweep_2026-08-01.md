# Reproduction sweep (2026-08-01, in progress)

Re-running every locally feasible experiment from its committed script and
fixed seeds, diffing against committed results. Strategy revised after a
45-agent parallel workflow hit the session limit with zero results: the
sweep now runs sequentially (inline in the main loop, or single small
agent batches when limits allow).

| dir | verdict | detail |
|---|---|---|
| 001_proxy-is-compression | REPRODUCED_FULL | all three identity checks pass; printed deviations match README magnitudes (1.3e-16 / 2.5e-16 / 1.1e-16), settling audit item U3 |
| 015_v2-density-law | REPRODUCED_FULL | rerun results.csv byte-identical to committed (280 instances); identities to 1e-10, cubic law <=3% |
| 014_argmax-robustness | REPRODUCED_FULL | rerun results.csv byte-identical (140 instances x 2 depths) |
| 041_sampled-plus-normalized | REPRODUCED_FULL | generated fresh from committed script during this loop (2026-07-31) |
| 042_sampled-N-mechanism | REPRODUCED_FULL | generated fresh during this loop |
| 043_synthetic-noise-control | REPRODUCED_FULL | generated fresh during this loop |
| 044_composed-recipe-n16 | REPRODUCED_FULL | generated fresh during this loop |
| 045_selection-aware-certificate | REPRODUCED_FULL | generated fresh during this loop; validity asserts held on every row |
| 016_cheap-prefix-frame | REPRODUCED_FULL | rerun results.csv byte-identical (210 jobs) |
| 033_statistics-hardening | ANALYSIS_VERIFIED | headline rhos/CIs recomputed from committed CSVs by the 2026-07-31 numbers audit |
| 034_ceiling-validation | ANALYSIS_VERIFIED | ceiling gaps + transfer regret recomputed from CSVs (audit; README refreshed to 0.006) |
| 035_conditioning-correction | ANALYSIS_VERIFIED | capture stats recomputed (min 0.9046, mean 0.97) |
| 036_regret-certificate | ANALYSIS_VERIFIED + addendum machine check rerun this loop (714 points) |
| 037_weighted-maxcut | ANALYSIS_VERIFIED | capture 92.9-99.6% recomputed from CSVs (audit) |
| 038_signed-regret-decomposition | ANALYSIS_VERIFIED | 557/560, 3.8-4.7x, margin ~40% recomputed from CSVs (audit) |
| 039_bias-corrected-objective | ANALYSIS_VERIFIED | rho range and 0/126 recomputed (audit) |
| 040_n16-stability | ANALYSIS_VERIFIED | 0.84/0.76 vs 0.25/0.47 convention split recomputed (audit) |
