# 020 — A two-point leakage certificate for regret, and a zero-cost normalization fix at depth

**Question:** The sharpest referee objection (journal_readiness gap 3): the
paper's error calculus bounds STATE error while its own E014/E033 show
parameter quality is an argmax-transfer matter, so "the theory bounds the
wrong error." Can the leakage calculus be connected rigorously to regret?

**Answer: Yes, by a two-point certificate (proved, and machine-verified on
all 280 instance-depth pairs), but it is honestly loose (median 10-11x), and
the looseness is itself the argmax-transfer diagnosis. Bonus finding:
computing the proxy objective from the NORMALIZED compressed state nearly
halves depth regret at zero cost (p=3: 0.080 -> 0.044, better on 134/140),
while p=1 mildly prefers the unnormalized convention (0.031 vs 0.047).**

## The certificate

With F the true normalized objective, F-hat the proxy objective from the
normalized compressed state chi = phi/||phi||, theta-star the true argmax
and theta-hat the proxy argmax over the same candidate set:

  regret <= eps(theta-star) + eps(theta-hat)

because F-hat(theta-star) <= F-hat(theta-hat) kills the middle term of the
three-term decomposition. And for any theta, inserting a = <C>_psi into
<psi|C|psi> - <chi|C|chi> = <psi-chi|(C-a)|psi> + <chi|(C-a)|psi-chi>
(valid since both states are normalized) and applying Cauchy-Schwarz:

  eps(theta) <= ||psi-chi|| (sigma_psi + sqrt(sigma_chi^2 + delta^2)) / c_opt,
  ||psi-chi|| <= ||psi-phi|| + (1-||phi||) <= 2 sum_l lambda_l(theta),

so the leakage calculus (Theorem 2) controls regret through leakage at just
two points, both of which sit in the low-leakage corner in practice.

## Measured tightness (results.csv, 140 instances x p in {1,3})

- The inequality held on every row (asserted in run.jl; a machine check of
  the proof).
- Median bound/regret ratio 10.5 (p=1) and 11.3 (p=3); range 6.4-46.
- Mean bound 0.46-0.48 AR against mean regret 0.044-0.047: the bound is
  informative only where accumulated leakage is a few percent (12/280 rows
  below 0.25 AR).
- Replacing the sigma factors by the trivial ||C|| = m makes it 10x worse
  (mean 4.7-5.5); replacing measured distance by the leakage sum costs a
  further 1.4-2.3x.

**Reading:** the calculus DOES bound the quantity of interest, and the
measured looseness quantifies exactly the paper's argmax-transfer thesis:
regret is small not because the proxy landscape is pointwise accurate but
because its errors at theta-star and theta-hat nearly cancel, which a
triangle inequality cannot see. The certificate is the honest bridge between
Theorem 2 and section 5.4, not a sharp predictor.

## The normalization finding

The certificate requires the normalized proxy objective, which changed the
selected grid point on 140/140 instances relative to the unnormalized
convention used by E002/E004/E014, and:

| depth | unnormalized regret | normalized regret | normalized better on |
|---|---|---|---|
| p=1 | 0.031 | 0.047 | 31/140 |
| p=3 | 0.080 | 0.044 | 134/140 |

Mechanism: along the p=3 ramp grid the compressed norm varies substantially,
so the unnormalized objective conflates "high predicted value" with "low
leakage" and is biased toward small-angle schedules; dividing by the tracked
norm (free, by Theorem 2's bookkeeping) removes the bias exactly where
leakage accumulates. At p=1 the bias is small and acts as a mild
regularizer, so the unnormalized convention keeps a small edge there.
Practical rule now in the paper: normalize the proxy objective at depth.

Reproduce: `JULIA_NUM_THREADS=auto julia --project research/experiments/036_regret-certificate/run.jl`
(~10 min); smoke: `E20_SMOKE=1 ...`.

## Addendum machine check (2026-07-31)

The unnormalized-objective extra term |a|(1-‖φ‖²)/c_opt (paper remark) was
audit-verified algebraically on 2026-07-31 and is now also machine-checked:
`verify_addendum.jl` asserts the Theorem-3 norm identity, the decomposition
identity, and the full unnormalized bound at 714 (instance, schedule) points
(7 families x n in {10,12} x 3 instances; p=1 grid and p=3 ramps), all to
1e-10 (bound strict; worst utilization 0.85).
