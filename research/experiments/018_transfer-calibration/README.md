# 018 — Transfer calibration: is the proxy better than just borrowing angles?

**Question.** How does the proxy's parameter-setting regret compare to the
standard cheap alternative — transferring angles grid-optimized on a small
(n=10) instance of the same family? (Paper-critique item T2.4: the one
calibration number a referee will predictably ask for.)

**Answer. Transfer wins every cell, usually by ~10×: pooled transfer regret
0.0072 (median 0.0026) vs 0.0606 for the exact-compression proxy; even the
worst transfer cell (dense ER, p=3, n=18: 0.042) beats the proxy's 0.09.**

Job 11677649 (65 s), 700 targets × 3 sources × 2 depths, all seed
regenerations asserted. Full table: [analysis.txt](analysis.txt)
(`julia analyze.jl`). Patterns worth noting:

- **Within-family parameter concentration is extremely strong at these
  sizes.** A single n=10 source's true-landscape argmax transfers to n=12–18
  targets at regret ≤ 0.005 on five of seven families (p=1); best-of-3
  sources drives 3-regular to ≈ 0.0003 (consistent with the fixed-angle
  literature) and WS(k=4;b=0.1) to ≈ 0.0001 (a near-deterministic family).
- **The one place transfer strains is dense ER at depth**: p=3 regret grows
  with n (0.008 at n=12 → 0.042 at n=18), the same regime where compression
  leaks most. Everywhere else transfer regret is flat or falling in n.
- **Implication for the paper (§5.4/§6):** at classically simulable sizes
  the proxy is dominated not only in speed (E013) but in regret by the
  cheapest standard method. This *supports* the argmax-transfer thesis
  (parameters concentrate, so a small instance's argmax is all you need)
  and sharpens the honest pitch: the proxy's unique regime is where no
  cheap source landscape exists — beyond-simulation sizes with the
  analytical N, and hardware settings where every evaluation costs shots.
  The understanding/error-calculus contribution is untouched by this.

## Method

For every target instance already measured in exps 002 (n=12,14; 30
instances/cell) and 010 (n=16,18; 20 instances/cell), across the seven
standard families: regenerate the instance exactly from the seed recorded in
those experiments' CSVs (edge count asserted against the CSV), then evaluate
real QAOA at angles grid-optimized on the *true* landscape of an n=10 source
instance of the same family — 3 independent sources per family — at p=1
(40×40 grid) and p=3 linear ramps (8⁴ endpoint grid), the same grids as
002/010. Transfer regret is measured against the target's own recorded grid
ceiling, so it is directly comparable, cell by cell, to the exact-compression
and sampled-N proxy regrets already in those CSVs. Because the source angles
live on the same schedule grids as the ceilings, transfer regret is
nonnegative by construction (asserted).

Output: `results.csv`, long format — one row per (target, source, p):
`family,n,instance,seed,source_rep,p,ar_transfer,ar_ceiling,regret_transfer`.

## Caveats

- Sources are n=10 with the *true* (statevector) landscape — the most
  favorable version of transfer (no surrogate error on the source side); this
  makes the comparison conservative *against* the proxy.
- Transfer is within-family with unscaled angles; rescaling rules
  (e.g. γ ∝ 1/√m) could improve the dense→dense transfers and are cited, not
  tested, here.
- p=3 transfers share the ramp-endpoint grid with the ceilings; finer
  schedules could shift both sides.

## Reproduce

```bash
cd research/experiments/018_transfer-calibration && sbatch run.sb
E18_SMOKE=1 julia --project research/experiments/018_transfer-calibration/run.jl  # smoke
```
