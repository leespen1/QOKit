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
