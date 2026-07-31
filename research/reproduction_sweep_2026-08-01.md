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
