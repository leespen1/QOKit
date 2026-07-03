# E014 — The fitted-shape paradox: is the weighted norm the right fit objective?

**Question.** Does Theorem 3's f_d(β)-weighted transfer error — not the
entrywise MSE that the G-RIPS shape fits minimized — control
parameter-setting regret?

**Answer.** *Pending — Slurm array job 11575420 running (submitted 2026-07-03).*

**Early signal from the n=10 smoke (single ER(0.5) instance, coarse grid):**
under the normalized ranking, the transfer-invisible `nulled` perturbation has
regret 0.0000 at every ε (entrywise MSE up to 0.23), while the aligned `lowd`
perturbation's regret grows 0.10 → 0.32 with ε — at matched entrywise MSE,
regret tracked the weighted error exactly. The raw Eq.-9 ranking was instead
destroyed by every perturbation (norm-inflation hijack, normw up to 4.7×10⁴ vs
0.978 for exact N — contractivity in action).

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
cd research/experiments/014_fitted-shape-paradox && sbatch run.sb   # full
E14_SMOKE=1 julia --project research/experiments/014_fitted-shape-paradox/run.jl  # smoke
```
