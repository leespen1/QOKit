# 018 — Is the grid ceiling a fair regret reference, and how does the proxy compare to parameter transfer?

**Question:** The journal-readiness review flagged (a) that regret is measured
against a grid ceiling never validated against continuous optimization
(especially the coarse 8^4 linear-ramp grid at p=3), and (b) that no external
baseline (published-style parameter transfer) is ever compared.

**Answer: (a) The ceilings are tight; the ramp restriction costs more than
the grid. (b) Parameter transfer from one brute-forced small cell BEATS the
exact-compression proxy on essentially every tested instance, decisively at
depth.**

## Results (results.csv, 140 instances: 7 families x n in {12,14} x 10, same
seeds as E002/E004/E014)

1. **p=1 ceiling is essentially exact.** Nelder-Mead refinement of the true
   objective from the 40x40 grid argmax raises the ceiling by mean 0.0003 AR
   (max 0.0007). No p=1 conclusion changes.
2. **p=3 ramp-grid ceiling is nearly tight within the ramp class.**
   Refinement over the 4 ramp endpoints gains mean 0.0023 (max 0.0049).
3. **The linear-ramp restriction itself costs more.** Releasing all 6 angles
   (NM from the refined ramp) gains a further 0.0067 mean, family-dependent:
   3-regular 0.019, ER(0.25) 0.010, dense families ~0.003; max 0.052 on one
   sparse instance. So p=3 "regret" is regret against the best coarse-grid
   linear ramp, understating regret against the unrestricted p=3 optimum by
   ~0.009 on average. The paper now states this.
4. **Transfer baseline (the sharp finding).** Averaging the true (grid)
   argmax angles of the ten ER(0.5), n=12 instances and applying them
   verbatim to all 140 instances (all families, both sizes):
   - p=1: transfer AR beats the proxy-set AR on 134/140 instances
     (mean advantage 0.018); pooled transfer regret ~0.014 vs proxy 0.032.
   - p=3: transfer beats the proxy on 140/140 (mean advantage 0.074);
     pooled transfer regret ~0.008 vs proxy 0.083.
   Parameter concentration on these random ensembles is strong enough that
   one solved small instance sets better parameters than instance-specific
   proxy compression, especially at depth (the ramp grid's coarse 8-point
   axes hurt the proxy argmax more than they hurt a transferred point).

## Reading

This does not contradict the paper's theory (the compression and its error
calculus are about state evolution, and E014/E017 already located parameter
quality in argmax transfer). It does sharpen the practical framing: on
concentrated homogeneous random ensembles at simulable sizes, per-instance
parameter setting (by any N) adds little over transfer from one solved
instance. The proxy's comparative advantage must lie where concentration
fails or no solved source exists (structured/heterogeneous instances,
untested here), or at sizes where nothing can be brute-forced and only the
analytical N survives. The paper's Discussion and abstract now say this.

## Method

Same instance generation as E014 (seed scheme SEED + 10000 fam + 100 n +
inst). True objective via `qaoa_expectation` statevector evaluation.
Box-clamped Nelder-Mead (400 iterations, ftol 1e-9) from the grid argmax;
6-angle refinement starts from the refined ramp. Transfer angles: mean of
the source cell's per-instance grid argmaxes, p=1 (gamma, beta) and p=3
ramp endpoints separately.

Reproduce: `JULIA_NUM_THREADS=auto julia --project research/experiments/018_ceiling-validation/run.jl`
(~15 min on 8 threads); smoke: `E18_SMOKE=1 julia --project .../run.jl`.
