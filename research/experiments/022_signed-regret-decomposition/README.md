# 022 — The exact signed anatomy of regret (why the certificate is 10x loose)

**Question:** E020's two-point certificate is honest but ~10x loose. Where
exactly does the slack come from? The paper's first reading ("errors at the
two points nearly cancel") deserved a measurement.

**Answer: regret = Delta-e - margin exactly, and three effects stack.
(1) The normalized proxy OVERPREDICTS the objective essentially everywhere
(e = F - F-hat < 0 at 277/280 measured argmax points), so the two signed
errors share a sign and partially cancel (median cancellation factor ~0.6,
correlation 0.86 at p=1). (2) A winner's curse inflates the error at the
proxy's own argmax: |e(theta-hat)| is 4-5x |e(theta*)| (0.10 vs 0.02-0.03),
because maximizing F-hat selects points of maximal overprediction. (3) The
proxy's own margin F-hat(theta-hat) - F-hat(theta*) absorbs ~0.03 of the
remaining ~0.077, leaving regret ~0.045. The certificate's Cauchy-Schwarz
step alone inflates each |e| by ~3.6x; dropping signs and margin costs the
remaining ~2.7x.**

## The identity (asserted to 1e-9 on all 280 instance-depth rows)

  regret = [e(theta*) - e(theta-hat)] - [F-hat(theta-hat) - F-hat(theta*)]
         =        Delta-e             -            margin,   margin >= 0.

## Measured (140 instances x p in {1,3}, normalized proxy objective)

| quantity | p=1 | p=3 |
|---|---|---|
| mean e(theta*) | -0.027 | -0.021 |
| mean e(theta-hat) | -0.104 | -0.097 |
| corr(e*, e-hat) across instances | 0.86 | 0.33 |
| same sign | 140/140 | 137/140 |
| mean Delta-e | 0.077 | 0.076 |
| mean margin | 0.030 | 0.033 |
| mean regret | 0.047 | 0.044 |
| median cancellation factor Delta-e/(|e*|+|e-hat|) | 0.58 | 0.61 |

## Reading

The refined story for the paper: regret is small for three stacked reasons,
none of which is pointwise landscape accuracy. The proxy's error field has a
systematic sign (uniform overprediction, which argmax comparison partially
sees through), the residual tilt is winner's-curse-concentrated at the
proxy's own argmax, and the proxy's internal margin absorbs a third of that
tilt. A sharp regret theory would need to model the selection effect, not
just the pointwise error, which is exactly why no pointwise norm predicts
regret (E012/E014).

Reproduce: `JULIA_NUM_THREADS=auto julia --project research/experiments/022_signed-regret-decomposition/run.jl` (~10 min).
