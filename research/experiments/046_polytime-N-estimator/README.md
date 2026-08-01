# E046 — Polynomial-time Monte Carlo estimator of N(c'; d, c)

## Question

Can a polynomial-in-n Monte Carlo estimator of the homogeneous distribution
N(c'; d, c) (plus estimated class counts) replace the exponential-cost
enumeration sampler in the sampled-N + normalized-objective recipe
(E041/E044) without losing parameter-setting quality? The paper states this
as an open question and confines the sampled-N recipe to n <= 30 because
`sampled_homogeneous_distribution` enumerates all 2^n bitstrings.

## Answer

**Yes.** The polynomial-time estimator (poly N + Wang-Landau counts +
normalized objective) matches the enumeration-sampled S=10 recipe within the
~0.01 AR success band, and at p=1 it is the best variant outright. Pooled
mean regret_norm over 7 families x n in {12,14} x 5 instances: p=1 poly
0.0112 vs samp10 0.0186 vs exact 0.0467; p=3 poly 0.0302 vs samp10 0.0249 vs
exact 0.0450. Paired per-instance (poly minus samp10): p=1 mean -0.0044
(n=12) / -0.0103 (n=14), poly better on 49/70; p=3 mean +0.0058 / +0.0046,
poly slightly worse but well within 0.01 and still roughly half the exact-N
regret. Wang-Landau missed 0 attained cost classes on all 70 validation
instances. Wall-clock crosses over as expected: the enumeration sampler is
faster at n <= 14 (~0.02s vs ~0.2s) but the poly estimator wins at n=20
(0.96s vs 1.70s mean, 3 ER(0.5) instances) and its cost is polynomial, so
the gap widens with n. This closes the paper's open question positively.

## Method

Estimator (`estimator.jl`, `polytime_homogeneous_distribution`) — uses only
the edge list, never the 2^n state space:

1. **Density of states**: Wang-Landau over cost classes c in 0..m (single-bit
   flip proposals, O(deg) incremental cost deltas; flatness 0.8 over visited
   classes, ln f halved from 1 to 1e-3, 2e6-step stage cap).
2. **Class members**: flat-histogram Metropolis production chain with
   stationary pi(x) proportional to 1/g(c(x)); thinned every n flips; up to
   S=10 distinct members per attained class (tiny classes yield fewer).
3. **Profile rows**: for member x and each distance d, K=200 random
   d-subsets of flip positions, cost evaluated from scratch (O(m)), scaled
   by C(n,d); distances with C(n,d) <= K enumerated exactly (Gosper's
   hack). Each row sums to C(n,d) exactly; each member profile to 2^n
   (asserted).
4. **Counts**: M_hat_c = 2^n softmax(ln g). The poly variant pairs its N
   with these estimated counts; exact/samp10 use exact counts.

Validation (`run.jl`): 7 families x n in {12,14} x 5 instances (seeds
SEED=20260611, instance_seed(fi,n,inst)=SEED+10000 fi+100 n+inst); variants
{exact, samp10, poly} x {raw, normalized} objective; p=1 grid (40x40) and
p=3 linear-ramp grid (8^4); regret vs the statevector grid ceiling
(harness cloned from E044). Timing: per-instance wall-clock of the
enumeration sampler vs the poly estimator, plus an n=20 scaling run
(3 ER(0.5) instances, estimator only).

## Result

Mean (max) regret_norm by cell, from `results.csv`:

| n  | p | exact           | samp10          | poly            |
|----|---|-----------------|-----------------|-----------------|
| 12 | 1 | 0.0435 (0.1195) | 0.0187 (0.0559) | 0.0142 (0.0528) |
| 12 | 3 | 0.0409 (0.1333) | 0.0224 (0.0640) | 0.0282 (0.0896) |
| 14 | 1 | 0.0499 (0.0873) | 0.0185 (0.0550) | 0.0082 (0.0317) |
| 14 | 3 | 0.0490 (0.1158) | 0.0275 (0.0773) | 0.0322 (0.0832) |

Per-family (pooled over n, p) poly is best or tied with samp10 in 5/7
families; ER(0.25) is poly's weakest (0.0337 vs samp10 0.0235). The raw
(unnormalized) objective is much worse for poly (pooled ~0.18 at p=3),
confirming E041's finding that the normalized objective is the essential
partner; estimated counts make normalization even more load-bearing.

Timing (`timing.csv`): enumeration sampler mean 0.019s / 0.016s / 1.70s at
n = 12 / 14 / 20; poly estimator 0.18s / 0.22s / 0.96s. Crossover between
n=14 and n=20; enumeration scales as O(S m 2^n), poly polynomially.

## Caveats

- Prototype-stage estimator lives in the experiment dir, not `src/`.
- Wang-Landau may miss cost classes of exponentially small measure; missed
  classes are recorded (`classes_missed` in timing.csv), and claiming an
  unattainable class is a hard assert. All 70 validation instances missed 0
  classes; at n=20 no exact reference was computed, so misses there are
  unchecked.
- The n=20 timing crossover rests on only 3 ER(0.5) instances and one
  hyperparameter setting (S=10, K=200, WL stage cap 2e6); constants, not
  the asymptotics, could shift with tuning.
- The p=3 paired comparison has poly slightly behind samp10 (mean +0.005,
  worst instance +0.061 at n=12); "match within ~0.01" holds on the mean,
  not per instance.
- Members are "distinct states from a thinned flat-histogram chain", not
  exact uniform draws without replacement within a class; residual
  autocorrelation is possible for tiny classes.
- Validation ceiling still needs the statevector, so validation sizes stay
  at n <= 14; only the estimator itself is exercised at n=20.

## Repro

```bash
JULIA_NUM_THREADS=auto julia --project research/experiments/046_polytime-N-estimator/run.jl   # ~30-45 min
E46_SMOKE=1 julia --project research/experiments/046_polytime-N-estimator/run.jl              # ~2 min
```
