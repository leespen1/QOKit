# 017 — Error directions: which norm on ΔN governs parameter-setting regret?

*(Formerly numbered E014; renumbered after the 2026-07-03 branch reconciliation
because the June-13 line owns experiments 012 and 014. Slurm job name and log
files keep the original e014/job-11575420 identifiers in sacct.)*

**Question.** Does Theorem 3's f_d(β)-weighted transfer error — not the
entrywise MSE that the G-RIPS shape fits minimized — control
parameter-setting regret?

**Answer. Yes — the weighted transfer error, not entrywise MSE, governs
regret (paradox resolved); but "fit in the weighted norm" does NOT rescue the
triangle/normal shapes (their representational error in the weighted norm is
irreducibly ≈ 1), and the normalized objective, while the right *measuring
instrument*, is a *worse parameter-setter* for the analytical model, whose
raw values turn out to be argmax-informative.**

Full numbers: [analysis.txt](analysis.txt) (regenerate with
`julia analyze.jl`). Slurm array 11575420, all 10 tasks COMPLETED, 150
instances × 23 variants = 3450 rows. Three headline results:

1. **Part A — the paradox dissolves.** At matched entrywise error
   (‖ΔN‖ = ε‖N‖), regret under the normalized ranking is governed by the
   weighted error: transfer-invisible profiles (`nulled`, `highd`) sit at the
   exact-N baseline regret (0.029–0.043 vs baseline 0.0415) *even at ε = 0.5*,
   while the aligned profile hits 0.13–0.17 already at ε = 0.01. Pooled
   Spearman(regret, ew_avg) = 0.74 vs Spearman(regret, mse_rel) = 0.38
   (per-cell: 0.66–0.84 vs 0.22–0.50). Correlations are capped by regret's
   grid quantization and early saturation of the aligned profile.
2. **Part B — entrywise MSE misranks; the weighted metric ranks.** The
   normal fit has *lower* entrywise error than the triangle fit (0.60 vs
   0.69 pooled) but *higher* regret in most instances — entrywise MSE orders
   the pair correctly in only 45/150 instances (worse than chance, the
   G-RIPS paradox exactly), while ew_avg orders it correctly in 105/150.
   Against PaperProxy the weighted metric is also the better predictor
   (104–118 of ~148 vs 88).
3. **Part C — the prescriptive fix FAILS (negative result).** Refitting in
   the weighted norm barely helps the triangle (77 better / 6 tie / 67 worse,
   mean Δregret +0.007) and hurts the normal (16/64/70, −0.016): even the
   weighted-norm-optimal shapes keep relative weighted error ≥ 0.92. The
   fitted-shape program fails on *representational* grounds — the shape
   families cannot express the transfer operator — not because G-RIPS chose
   the wrong fit objective. This independently supports the paper's sampled-N
   recipe (exp 010), which needs no shape class at all.

**The raw-objective sideshow (norm inflation).** Thm-1 contractivity held on
every instance (exact-N norm weight ≤ 0.984 ≤ 1), while *every* coherent
perturbation — even ε = 0.01, transfer-invisible — inflated some grid point's
norm enough to hijack the raw Eq.-9 argmax (pooled raw regret ≈ 0.16 for all
16 perturbation variants regardless of shape or scale; median norm weights up
to 3.8×10⁵ for fitted shapes). But normalization is NOT a better
parameter-setting recipe: it slightly hurts exact N (raw 0.0315 vs nrm
0.0415, raw better in 99/150) and badly hurts the analytical PaperProxy
(raw 0.051/0.045 → nrm 0.16; raw better in ~132/150). The exp-004–006
finding must be refined: the analytical proxy's values carry no usable
*absolute* information, but its raw (unnormalized) landscape is
argmax-informative — its norm inflation is *correlated with parameter
quality* off the dense-ER pathology, so dividing it out discards signal.
The practical recommendation is unchanged (sampled N + empirical P + raw
objective); the normalized objective is the right instrument for *measuring*
model-error effects (Parts A–C above), not for setting parameters.

## Why this experiment

The old research log claimed (hypothesis H2, never re-verified): triangle and
normal proxies fitted to lower entrywise MSE against the empirical N(c';d,c)
gave *worse* parameter setting. Theorem 3 offers an explanation: the proxy
update only sees N through the contraction T(c',c;β) = Σ_d f_d(β) N(c';d,c)
with f_d(β) = cos(β)^(n−d)(−i sin β)^d, so error components along d-profiles
where f is small are invisible, and entrywise MSE weights them all equally.
The parameter-relevant model error is the transfer-matrix (Frobenius) error

E_w(β)² = Σ_{c',c} |Σ_d f_d(β) ΔN(c';d,c)|²

(γ phases are unimodular per column and drop out). Global scale is a gauge
freedom of the proxy argmax, so all metrics and fit objectives optimize a
global scale factor in closed form first.

## Method

All parts run at p=1 on the standard 40×40 (γ,β) ∈ [0,π]×[0,π/2] grid with
true-QAOA grid ceilings, on n ∈ {12, 14}, 15 instances × 5 families
(ER(0.5), ER(0.25), BA(k=4), WS(k=4;b=0.1), 3-regular), seeds fixed —
mirroring experiments 002/010. Regret = AR_ceiling − AR at the
proxy-chosen grid point, evaluated by real statevector QAOA.

**Two argmax rankings are recorded for every variant.** The smoke run
exposed a degeneracy in the paper's raw Eq.-9 objective
(⟨C⟩ = Σ 2ⁿP(c)|Q(c)|²c, used by exps 001–010): *any* coherent model
perturbation, however small (ε = 0.01, transfer-invisible), plants a
norm-inflated beacon somewhere on the grid (predicted ⟨C⟩ growing like ε²,
up to ~5×10⁵ on a 21-edge graph) that hijacks the raw argmax — the same
pathology exps 004–006 diagnosed for the analytical proxy on dense ER. So
each variant also gets the **normalized** ranking: the same quantity divided
by the state weight Σ 2ⁿP(c)|Q(c)|² (⟨C⟩ of the normalized compressed
state; gauge-invariant, immune to pure norm inflation). The normalized
ranking is the instrument for Parts A–C (it isolates landscape-shape error);
raw-vs-normalized disagreement is itself data — exps 005/006 showed
norm *vetoes* fail, and this tests whether norm *division* is the recipe
that works.

- **Part A (controlled perturbations).** Perturb the exact empirical N by
  ΔN[c',d,c] = u_d·‖N[c',:,c]‖ scaled to ‖ΔN‖_F = ε‖N‖_F,
  ε ∈ {0.01, 0.05, 0.2, 0.5}, with four d-profiles u: `lowd` (aligned with
  Re f(β=0.3) — maximally visible), `highd` (mass at d≈n/2, invisible at
  small β), `nulled` (an exponentially
  low-d-concentrated bump with span{Re f, Im f} at β ∈ {0.15, 0.3, 0.45}
  projected out — invisible at those βs by construction), `random`. At fixed ε the entrywise
  MSE is identical by construction; prediction: regret tracks E_w, spanning
  orders of magnitude across profiles.
- **Part B (fits as G-RIPS did).** TriangleProxy (4 params, via
  `IntuitiveTriangleProxy`) and NormalProxy (3 params) fit to the exact
  empirical N by gauge-fixed entrywise error, using seeded random search
  (1500 uniform samples + 4 shrinking-box rounds × 400). Plus the analytical
  PaperProxy with effective edge probability 2m/(n(n−1)), evaluated with both
  empirical P (`paper:empP`) and its own binomial P (`paper:binP`). Question:
  when entrywise MSE misranks models vs regret, does E_w rank correctly?
- **Part C (prescriptive fix).** Fit the same two shapes minimizing E_w
  averaged over 8 reference βs spanning the grid instead. Prediction:
  weighted-norm fits set parameters at least as well as entrywise fits
  despite much larger entrywise MSE — turning the diagnosis into a recipe
  ("fit in the norm the dynamics sees").

Output: `results_task<ID>.csv`, long format, one row per (instance, variant):
gauge-fixed `mse_rel`, `ew_avg` (over the grid's βs), `ew_star` (at the
exact-N proxy's chosen β), real-QAOA `ar`, `regret`, argmax displacement
(`dgamma`, `dbeta`), the proxy's raw predicted ⟨C⟩ (`pred_exp`, norm-inflation
diagnostic), and fit parameters.

## Relation to experiment 012 (the June-13 independent test)

Exp 012 tested the same H2 question with the actual G-RIPS fitting
conventions (slice-normalized MSE, `fit_proxy_to_real` optimizer settings,
fit-trajectory snapshots) and found that **no scalar mismatch norm — not
entrywise MSE, not one-layer amplitude error, not landscape correlation —
predicts regret across its 10-model zoo; only argmax displacement does
(ρ≈0.7)**. This experiment's Part A explains why both results are right: the
weighted transfer error *is* the governing quantity when it has dynamic range
(controlled directions span 4 orders of magnitude in it at fixed MSE), but
every realistic shape fit *saturates* it (≥ 0.92 here even after
weighted-norm fitting; 012's zoo likewise all-bad), so within a model class
the norm carries no discriminating signal and the downstream argmax
displacement becomes the operative predictor. Two design differences also
matter: 012's amplitude norms were scale-hostage (raw PaperProxy error ~1e5
from slice-sum conventions) whereas the transfer-matrix metric here is
gauge-fixed — global scale optimized out, which is exactly the invariance of
the argmax; and 012 measured state-level error while E_w is operator-level.
Joint conclusion for §6: *direction of model error relative to the dynamics,
not its size, is causal; shape families err almost entirely in visible
directions; argmax displacement is the only within-class discriminator; and
no fit objective rescues the shapes.*

## Figures

*(added after analysis)*

## Caveats

- Perturbed N arrays may have (small) negative entries; the proxy iteration
  is linear and does not require positivity. This is a synthetic stress test
  of the *metric*, not a physical model class.
- The `nulled` profile is exactly invisible only at its three reference βs;
  `ew_avg` over the full grid is small but nonzero.
- Fits use the exact empirical N of the *same instance* as the target, i.e.
  the best case for shape models — model-class error only, no
  instance-transfer error.
- `pred_exp` is reported without renormalization; only its argmax matters to
  parameter setting (global scale is gauge).

## Reproduce

```bash
cd research/experiments/017_error-directions && sbatch run.sb   # full
E17_SMOKE=1 julia --project research/experiments/017_error-directions/run.jl  # smoke
```
