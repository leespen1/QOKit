# E046 — Polynomial-time Monte Carlo estimator of N(c'; d, c)

## Question

Can a polynomial-in-n Monte Carlo estimator of the homogeneous distribution
N(c'; d, c) (plus estimated class counts) replace the exponential-cost
enumeration sampler in the sampled-N + normalized-objective recipe
(E041/E044) without losing parameter-setting quality? The paper states this
as an open question and confines the sampled-N recipe to n <= 30 because
`sampled_homogeneous_distribution` enumerates all 2^n bitstrings.

## Answer

**PENDING (full run in progress).** Smoke (n=10, 2 instances/family):
poly+normalized mean regret 0.0145 (p=1) / 0.0171 (p=3) versus samp10+normalized
0.0179 / 0.0123 and exact+normalized 0.0264 / 0.0214 — well within the 0.01
match band of the enumeration recipe.

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

PENDING. Smoke results in `results_smoke.csv` / `timing_smoke.csv`.

## Caveats

- Prototype-stage estimator lives in the experiment dir, not `src/`.
- Wang-Landau may miss cost classes of exponentially small measure; missed
  classes are recorded (`classes_missed` in timing.csv), and claiming an
  unattainable class is a hard assert. Smoke missed 0 classes everywhere.
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
