# Research status

*Last updated: 2026-07-05. This is the one page Spencer needs to read.
Everything here links to a reproducible experiment or a committed document.
Plain-language walkthrough of all results: [explainer.md](explainer.md).
Running journal of the autonomous work: [journal.md](journal.md).*

## ⚠ Branch reconciliation (2026-07-03)

Two autonomous work lines had diverged from commit `8680853` without knowing
about each other: an **overnight session of 2026-06-13** (experiments 012,
014–016, the literature pass, the paper critique — pushed to
`origin/ClaudeResearch` but never pulled here) and the **July sessions** on
this clone (E013 write-up, explainer/journal, and a fitted-shape experiment
that unknowingly re-tested exp 012's question with a different design). They
were merged today. Consequences:

- The July fitted-shape experiment is **renumbered 017** (the June line
  already owned 012 and 014); its design and exp 012's turn out to be
  *complementary, not redundant* — see the joint entry below.
- The June-13 session's **paper-repo edits were committed only on an
  unpushed clone and are presumed lost** (§5.3 ranking figure + regret
  table, three draft fixes, two number-audit corrections, five citations).
  All of their *content* is preserved in committed QOKit files
  (`paper_critique.md`, `literature_pass.md`, `proposed_paper_additions.md`,
  and the figure generators in exps 004/011/012) and is being re-applied to
  the paper.
- `explainer.md` (results walkthrough, July line) and `paper_explainer.md`
  (paper companion, June line) both exist; they serve different purposes but
  overlap — consolidation is queued as a Spencer decision.

## Overnight digest (2026-06-13, from the merged line)

1. **`research/paper_explainer.md`** — plain-language companion to the paper:
   the proxy-as-compression reframe, the two error axes, the theorems as a
   story.
2. **`research/paper_critique.md`** — adversarial review + full number audit
   (every quoted statistic checked; two errors found). *[Paper-side fixes
   lost with the unpushed clone — being re-applied.]*
3. **§5.3 ranking figure + regret table** — generator committed (exp 004);
   the figure shows the leakage→regret ranking honestly (ρ=0.96/0.86 at p=1,
   ≈0 at p=3). *[Inclusion in the tex lost — being re-applied.]*
4. **New result E015** — the V₂ black box behind the density law opened: two
   exact lemmas + a triangle-conditioning mechanism; paper-ready LaTeX in
   `research/v2_density_law.md`.
5. **New result E014 (argmax-robustness)** — depth regret is argmax-transfer,
   decoupled from fidelity at p=3; resolves the depth framing question.
6. **New result E016** — cheap-prefix low-rank frames can't find the moving
   trajectory subspace (71→88° rotation); the Discussion's chicken-and-egg
   holds. A measured obstacle, not a method.
7. Five decisions for Spencer surfaced (see bottom).

**Goal:** QCE 2027 contributed paper, *"When and why does the homogeneous proxy work?
QAOA parameter setting as subspace compression."* Complete draft by **Aug 15, 2026**;
arXiv late August; submit April 2027. Full plan: see `research/decisions.md` and
`research/program.md`.

## Where we are

**Phases 1–4 effectively complete.** The paper
(`papers/OverleafPaper`, branch `ClaudeResearch`, head `7cdf80c`) compiles
clean at 10 pp with the theory (V₂ lemmas included), the full anatomy of
experiments 001–028, eight figures (ranking, argmax-transfer,
direction-vs-size, weighted binning, α-crossover with pre-registered
points), a glossary, a Conclusion, and 36 verified references from two
independent literature passes. Four adversarial review passes have been
applied (June 13 critique; July 4 referee; July 5 coherence; July 5 final
verdict: **weak accept, mechanical fixes only — all applied**). The
practitioner map is measured in every dimension: transfer wins
family-concentrated regimes (unweighted, mild weights, all sizes/depths
tested); the per-instance binned proxy wins high-dispersion weights at low
depth (boundary at order-unity CV², pre-registered and confirmed;
scale-checked to n=26); the deep high-dispersion corner is hard for every
method and marked open; Max-3-XOR replicates the map. Remaining: author
list + acknowledgments (parked per Spencer), the reviewer's optional
length-cut menu, and Spencer's read. Overleaf shows only `main`; merge the
paper repo's ClaudeResearch branch when ready to edit there.

## What we know (established results only)

- **The homogeneous proxy is exactly an orthogonal subspace compression of QAOA**
  (not an approximation): running the proxy with the same-instance empirical
  N(c';d,c) equals evolve-one-layer-then-project-onto-cost-classes, to machine
  precision. → [experiments/001_proxy-is-compression](experiments/001_proxy-is-compression/README.md)
  Consequence: total proxy error splits cleanly into *compression error*
  (leakage out of the cost-class subspace, model-independent) + *model error*
  (using analytical/fitted N instead of the instance's own).
- **Parameter setting is well-posed, and the compression itself is nearly lossless
  at p=1 on every graph family** (ER, BA, WS, 3-regular; n=12–14, 420 instances):
  value-added over the balanced-partition baseline is positive for every single
  instance, and the exact-compression proxy lands within ~0.03 AR of the true grid
  ceiling on all families. Regret grows to ~0.06–0.09 at p=3 ramps.
  → [experiments/002_baselines-and-headroom](experiments/002_baselines-and-headroom/README.md)
  Two consequences: (i) the old log's "proxies fail on non-ER graphs" must be a
  *model-error* story (bad analytical N), not a compression-error story, at least
  at p=1 — sharpens what E1.3 must separate; (ii) depth, not graph family, is
  where compression error first bites.
- **The Theorem-2 bound is tight enough to use, and leakage predicts fidelity.**
  At p=20 ramps (840 runs, n=12–16): ‖ψ−φ‖ ≤ Σλ holds with median slack ≈ 4×
  (never >10×), Σλ is a near-functional predictor of the actual distance, and it
  ranks families by fidelity at Spearman ρ ≈ 0.96–1.0 (except the decayed
  n=16/large-angle cell). → [experiments/003_leakage-vs-overlap](experiments/003_leakage-vs-overlap/README.md)
- **Surprise: density, not "ER-ness," drives compression error** — family-mean Σλ
  tracks edge count at Pearson 0.97; ER(0.5) is the *worst*-compressing family
  tested, 3-regular the best (at fixed unscaled angles; see exp 003 caveats).
- **GATE (exp 004): leakage ranks families by parameter-setting regret at p=1**
  (Spearman ρ = 0.86–0.96, criterion ≳0.8 met); at p=3 along proxy-chosen ramps
  both quantities barely vary across families (ρ ≈ 0 — no signal, not a
  contradiction). → [experiments/004_gate-leakage-vs-regret](experiments/004_gate-leakage-vs-regret/README.md)
- **H-ER-specific (old log) is contradicted**: analytical PaperProxy with effective
  edge probability gets regret ≈ 0.01–0.02 at p=1 on BA/WS/3-regular/sparse-ER —
  *better* than exact compression (the smooth model regularizes parameter choice).
  Its one failure is dense ER(0.5), where its argmax lands on an *unphysical
  norm-inflated artifact* (predicted ⟨C⟩ = 93 on a 38-edge graph). Since exact
  compression is contractive (Thm 1), norm inflation certifies model error for
  free; a norm sanity filter repairs the diagnosed instance
  (`004/diagnose_paper_artifact.jl`); exps 005/006 quantify why value-based
  filters nonetheless fail as a recipe.
- **Timing (exp 013): the exact N is the pipeline's wall, and brute force is
  free at simulable sizes.** On an A100, exact N costs 21.8 s at n=20 (clean
  O(4ⁿ) growth) vs 0.48 s for the brute-force p=1 ceiling it competes with;
  sampled N (S=10) is 50× cheaper (0.44 s). GPU statevector passes are
  kernel-launch-bound (~0.2–0.3 ms, flat in n ≤ 20), so exhaustive grid search
  matches the proxy pipeline end to end even at p=20; the pipeline's cost
  advantage opens around n ≈ 22–24 (extrapolated). Running n=20 exact N
  required grid-striding the GPU homodist kernel (commit `ce337950`). The
  paper prints the measured table in §5.3.
  → [experiments/013_timing-benchmark](experiments/013_timing-benchmark/README.md)
- **V₂ (the cubic-leakage variance functional) has an exact anatomy (exp 015).**
  Two machine-verified lemmas (1e-10, 280 instances): (L1) V₂ is the within-cost-class
  variance of T(y)=Σ_i(Σ_{j∼i}s_j)², and (L2) its *unconditional* variance is exactly
  2·Σ_{j≠k}A_{jk}² (squared codegrees), → 8p²m² for ER. The cubic law λ₁=(βγ²/8)√V₂
  holds to <0.3% everywhere. The density law √V₂∝m is an approximate cancellation:
  codegree variance grows super-linearly with density, conditioning on cost removes a
  triangle-driven fraction (ρ≈0.37 dense ER → 0.80 sparse 3-regular), and the two
  nearly cancel. Advances the open V₂ problem to one quantity (the conditioning
  correction). → [experiments/015_v2-density-law](experiments/015_v2-density-law/README.md);
  derivation + paper-ready LaTeX in [research/v2_density_law.md](v2_density_law.md)
- **Parameter-setting regret is argmax-transfer, decoupled from fidelity at depth
  (exp 014).** Pooled over 140 instances, regret tracks argmax displacement at both
  depths (ρ=0.76 at p=1, 0.63 at p=3) but the fidelity deficit only at p=1
  (ρ=0.39→−0.02 at p=3). Landscape flat-peak robustness does not predict regret
  (ρ≈−0.18/+0.10). The leakage calculus bounds the *state* error; parameter
  *quality* is a separate parameter-space matter — direct evidence for §5.4 and the
  mechanism behind "a worse-fidelity model picks better parameters."
  → [experiments/014_argmax-robustness](experiments/014_argmax-robustness/README.md)
- **An instance-adapted low-rank frame is not cheaply discoverable (exp 016).** The
  p=20 ramp trajectory is ~4-dim (oracle PCA captures ~0.99, confirming exp 009),
  so a good 4-dim frame would beat the (m+1)-dim cost-class frame — but a frame
  built from a cheap 5-layer prefix is near-orthogonal to the true late subspace
  (principal angle 71→88° growing with ramp) and captures *less* than the zero-cost
  cost-class frame. The moving subspace rotates over depth; the chicken-and-egg of
  the Discussion's "most interesting open question" holds. A measured obstacle, not
  a method. → [experiments/016_cheap-prefix-frame](experiments/016_cheap-prefix-frame/README.md)
- **The fitted-shape paradox is fully dissected — two independent designs agree
  and interlock (exps 012 + 017).** exp 012 (June 13, the G-RIPS fitting
  conventions): entrywise-MSE fitting is harmful (Triangle: regret 0.054 → 0.211,
  worse on 136/140) or inert (Normal: 10× better MSE, argmax never moves), and
  across the fitted-model zoo *no scalar mismatch norm* predicts regret — only
  argmax displacement does (ρ≈0.7). exp 017 (July 3, controlled error
  *directions* at matched entrywise MSE, gauge-fixed metrics): Theorem 3's
  f_d(β)-weighted transfer error *does* govern regret when it has dynamic range
  (pooled ρ=0.74 vs 0.38 for MSE; transfer-invisible perturbations at 50%
  entrywise error carry zero excess regret) — but every realistic shape fit
  saturates that norm (≥0.92 *even when fitted in it*), which is exactly why
  012 saw no discriminating scalar across models. Joint story for §6: error
  *direction* relative to the dynamics is the causal quantity; shape families
  err almost entirely in visible directions, so no fit objective can save them
  (017 Part C), and within a model class only argmax displacement
  discriminates (012, 014). Bonus: *any* coherent model error hijacks the raw
  Eq.-9 argmax via norm inflation (predicted ⟨C⟩ up to 5×10⁵; Thm-1
  contractivity held on all 150 instances, exact-N weight ≤ 0.984), yet
  dividing the norm out *hurts* the analytical proxy (regret 0.05 → 0.16) —
  its raw landscape is argmax-informative, sharpening 012's "amplitude norms
  are scale-hostage" and the 004–006 "values are noise" claim.
  → [experiments/012_fitted-shape-paradox](experiments/012_fitted-shape-paradox/README.md),
  [experiments/017_error-directions](experiments/017_error-directions/README.md)

## Working hypotheses (NOT established — from the deleted research log or intuition)

- H1: The proxy's usefulness is governed by leakage out of the cost-class subspace,
  which concentrates only for ER-like graphs. (The paper's central claim.)
- ~~H2: Fitted proxy shapes (triangle/Gaussian) with lower entrywise MSE against the
  empirical N(c';d,c) gave *worse* parameter setting.~~ **Confirmed twice and
  explained** (exp 012: fitting harmful or inert under the G-RIPS conventions;
  exp 017: entrywise MSE is an *anti-predictor* — it orders triangle-vs-normal
  against regret in 105/150 instances — while Theorem 3's weighted transfer
  error is the quantity that governs, and shape families saturate it).
- H3: Sampling 5–10 bitstrings per cost class suffices to estimate the quantities that
  matter. (Old log; needs re-verification in the leakage metric, not entrywise N.)
- ~~H4: On ER(0.5), random balanced partitions reach ~75–85% approximation ratio.~~
  **Confirmed** (exp 002: 0.75–0.76 mean on ER(0.5) and BA(k=4) at n=12–14); folded
  into metric design as "value-added" = AR_proxy − AR_baseline.
- H2/H-ER-specific is now under pressure: exp 002 shows the *compression* is
  near-lossless on non-ER families at p=1, so any non-ER failure must come from the
  analytical/fitted N (model error). E1.3 tests this directly.

## Open questions / next experiments

- ~~E0.1: is the proxy numerically identical to the compressed evolution?~~ **Done, yes** (exp 001).
- ~~E1.1: is there headroom above baselines?~~ **Done, yes, on every family** (exp 002).
- ~~E1.2: does Σλ track the overlap deficit?~~ **Done, yes — bound slack ~4×,
  near-functional predictor** (exp 003).
- ~~E1.3 gate~~ **Run: PASS at p=1 (ρ = 0.86–0.96); p=3 inconclusive** (exp 004).
- ~~E005: norm filter as universal recipe~~ **Done: NEGATIVE — rescues only dense
  ER(0.5) (regret 0.159 → 0.077) and is catastrophic elsewhere; the analytical
  model's calibration breaks long before its argmax moves.** The argmax-robustness
  of the analytical proxy is itself a finding for §6.
  → [experiments/005_norm-filtered-paper-proxy](experiments/005_norm-filtered-paper-proxy/README.md)
- ~~E006: physicality cap~~ **Done: NEGATIVE — the model inflates predictions even
  at its (correct) peak on sparse families, so any value-based veto rejects the
  answer; on dense ER(0.5) the spurious region extends below the physical cap.**
  Filter arc (004→005→006) closed: the analytical proxy's values are pure noise in
  absolute terms; its argmax location is the only usable signal — excellent off
  dense graphs, corrupted on dense ER(0.5) by a spurious large-β peak. (Refined by
  exp 017: the raw value *landscape* is argmax-informative; see What-we-know.)
  → [experiments/006_physicality-filter](experiments/006_physicality-filter/README.md)
- ~~E007 = E2.1: leakage anatomy~~ **Done — three results:** (i) small-angle law
  λ ≈ const·β·γ²·m, with the O(βγ) term killed by an exact MaxCut identity
  (Σ_i c(x⊕e_i) = (n−4)c(x) + 2m — now a committed test) → quantitative
  explanation of the proxy's small-γ accuracy, lemma for §3/§4; (ii) density law
  sharpened: λ/(βγ²m) constant to ±40% across families; (iii) compression leakage
  does NOT flag the analytical proxy's spurious peak (38th percentile) —
  compression and model error are empirically independent axes, a central claim.
  → [experiments/007_leakage-anatomy](experiments/007_leakage-anatomy/README.md)
- ~~E2.2 sampled-leakage predictor~~ **Done: S=5 estimates η_F to median 3.2%
  error (uniform across families/angles) — far below the between-family
  differences; Phase-3 scale-up enabler confirmed. Bonus: Theorem 3's identity
  verified to 4.5×10⁻¹⁵ over 210 full-enumeration rows.**
  → [experiments/008_sampled-leakage-predictor](experiments/008_sampled-leakage-predictor/README.md)
- ~~E2.3 trajectory-PCA~~ **Done: "wrong subspace," decisively — p=20 ramp
  trajectories have effective dimension 2–4 (99% energy), so a 2–3-dim PCA
  subspace matches the entire (m+1)-dim cost-class subspace; proxy degradation
  is mis-aim of the fixed cost-class frame, not state complexity.**
  → [experiments/009_trajectory-pca](experiments/009_trajectory-pca/README.md)
- **E3.1 (scale-up, n=16–18): the structure survives scale unchanged** — regret
  flat in n (p=1: 0.03–0.05; p=3: 0.08–0.11), value-added positive everywhere at
  family level, and **the sampled-N proxy (S=10) matches exact-N parameter
  choices within ~0.01 AR with slightly LOWER regret** — the practical recipe
  the filter arc failed to find (sampled N + empirical P). Family regret
  differences compress at scale, so ranking claims must be scoped.
  → [experiments/010_scaleup-ranking](experiments/010_scaleup-ranking/README.md)
- **E3.2 (depth, p=30, n=16–20): leakage accumulates linearly in depth,
  quantitatively as Theorem 3 predicts** (ratio 1.505 vs 1.5 at small ramps),
  sublinearly in m along deep trajectories, with p=30 small-ramp overlaps of
  0.71–0.81 even at n=20.
  → [experiments/011_depth-scaling](experiments/011_depth-scaling/README.md)
- ~~E2.4 fitted-shape paradox~~ **Done: H2 confirmed strongly, and the planned
  explanation is refuted in its simple form.** Entrywise-MSE fitting of shapes is
  harmful (Triangle: mean p=1 regret 0.054 → 0.211, worse on 136/140) or inert
  (Normal: 10× better MSE, argmax identical to the unfitted default on 140/140).
  No scalar mismatch norm tested predicts regret (entrywise MSE ρ≈−0.1,
  one-layer amplitude error ρ≈−0.2, landscape Pearson ρ≈+0.1); argmax
  displacement does (ρ≈0.7). Raw PaperProxy has amplitude error ~1e5 (slice sums
  up to 9e8) yet near-best regret — amplitude norms are hostage to scale
  conventions the argmax ignores. Per-slice calibration to 2^n is argmax-neutral
  at p=1 (0/140 changes for shapes). Extends the 004→005→006 arc: for *every*
  model class, argmax location is the only usable signal; §6 recommendation is
  "don't fit shapes by MSE — use defaults or analytical N."
  → [experiments/012_fitted-shape-paradox](experiments/012_fitted-shape-paradox/README.md)
- ~~E029 corner refinement~~ **Done (2026-07-05, array 11773952): E028's
  hard corner is a GRID ARTIFACT — full-angle compass refinement from any
  start (proxy/transfer/universal/grid-ceiling) converges to one optimum
  (residual 0.000–0.006); even the 8⁴ ceilings sat 0.02–0.03 low. The
  deep high-dispersion landscape is sharp, not trap-riddled: cheap local
  refinement erases the corner wherever true evaluations are affordable.**
  → [experiments/029_corner-refinement](experiments/029_corner-refinement/README.md)
- ~~E028 dispersion at depth~~ **Done (2026-07-05, array 11770562): the
  niche narrows at depth but survives under lognormal (proxy 0.05–0.11
  vs every single source 0.07–0.29; pooled sources recover on sparse);
  under Pareto at p=10–20 all methods degrade toward parity — the deep
  high-dispersion corner is hard for everything tested. Map complete in
  all dimensions measured.**
  → [experiments/028_dispersion-at-depth](experiments/028_dispersion-at-depth/README.md)
- ~~E027 Max-3-XOR replication~~ **Done (2026-07-05, array 11757023): the
  map generalizes — compression regret at MaxCut magnitude, sampled ≥
  exact (4th occurrence), transfer/universal near-perfect on unweighted
  3-XOR, λ ∝ clause density; only the cubic cancellation is
  MaxCut-specific, as proved. Paper's Limits updated.**
  → [experiments/027_max3xor](experiments/027_max3xor/README.md)
- ~~E026 pre-registered prediction test~~ **Done (2026-07-05, array
  11756459): the dispersion theory's called shot lands on dense ER — flip
  at α≈2.4 (CV²≈1.0), inside the window written down before the run;
  3-regular keeps transfer at α≥2.2 (its flip is below α=2), so the
  order-unity threshold has a family-dependent constant.**
  → [experiments/026_predicted-flip](experiments/026_predicted-flip/README.md)
- ~~E025 lognormal separation~~ **Done (2026-07-05, array 11755955): the
  mechanism is dispersion, not variance-finiteness — lognormal weights
  (finite variance at every σ) reproduce the crossover at CV² ≈ 2–8;
  transfer collapses to 0.16–0.35 by σ=2.5 while the proxy improves with
  dispersion (to 0.004). Paper's boundary claim upgraded to "weight
  dispersion of order unity; infinite variance sufficient, not
  necessary" (`8337edd`).**
  → [experiments/025_lognormal-separation](experiments/025_lognormal-separation/README.md)
- ~~E024 heavy-tail at scale~~ **Done (2026-07-05, array 11754145): the
  α-boundary survives n=22/26 with the sampled-only proxy — α=1.5: proxy
  wins all cells (0.015–0.025 vs 0.030–0.130); α=3: transfer wins all.
  The positive niche is scale-robust and delivered by the practical
  variant.**
  → [experiments/024_heavy-tail-at-scale](experiments/024_heavy-tail-at-scale/README.md)
- ~~E023 tail-exponent sweep~~ **Done (2026-07-05, array 11751921): the
  transfer↔proxy crossover sits at the infinite-variance boundary α ≈ 2 —
  the proxy's regret is nearly α-independent while transfer sweeps through
  it (wins α>2, ties α=2 on dense ER, loses α<2). Principled, not tuned;
  in the paper (`5302abc`).**
  → [experiments/023_tail-exponent-sweep](experiments/023_tail-exponent-sweep/README.md)
- ~~E022 heavy-tail crossover~~ **Done (2026-07-05, array 11751040): the
  proxy's measured niche exists. Pareto(1.5) weights defeat mean-weight-
  rescaled transfer (the mean is tail-dominated) and the binned
  per-instance proxy wins every p=1 cell by 2–7× (0.012–0.022 vs
  0.031–0.111), staying stable at p=3 while transfer sources become a
  lottery (0.021–0.140). Exp(1) stays transfer territory. Practitioner
  map complete; paper's abstract/contribution 5/practical/Conclusion all
  updated (`7613c37`).**
  → [experiments/022_heavy-tail-weights](experiments/022_heavy-tail-weights/README.md)
- ~~E021 weighted transfer rematch~~ **Done (2026-07-05, array 11748691):
  no crossover — weights erode transfer 3–10× (and collapse the
  zero-knowledge universal schedule on dense ER at p=3) but same-family
  weighted sources still beat the binned proxies in every cell
  (0.003–0.022 vs 0.033–0.111). The proxy's per-instance steadiness is
  real; i.i.d. U[0,1] weights aren't heterogeneous enough to flip the
  ranking. In the paper alongside E020.**
  → [experiments/021_weighted-transfer](experiments/021_weighted-transfer/README.md)
- ~~E020 binned weighted MaxCut~~ **Done (2026-07-04, jobs 11744099 +
  11745830): the top referee risk became a result. Quantile-binned classes
  obey λ² ≈ λ_struct² + O(1/K); the K-dim binned proxy recovers
  integer-cost regret by K≈64–128 (3-regular n=16 p=3: 0.103 vs 0.105
  unweighted) at O(K²n) cost independent of m. In the paper (abstract,
  practitioner section, Limits rescoped).**
  → [experiments/020_binned-weighted](experiments/020_binned-weighted/README.md)
- ~~E019 claimed-regimes test~~ **Done (2026-07-04, GPU array 11740470):
  both untested regimes go to transfer. At p=10/20 ramps every proxy
  variant collapses (exact 0.07–0.14, sampled 0.11–0.21, analytical
  0.18–0.37) while transfer reaches 0.001–0.03 — even the zero-model-error
  compression misplaces the argmax by ~0.1 AR while its fidelity is still
  0.7–0.8. At n=22/26 (beyond exact N) transfer wins every cell including
  dense-ER p=3 (0.061 vs sampled 0.078+). Positives: analytical argmax
  improves with n on sparse regular graphs; its binomial P is load-bearing
  (empirical P: 0.05 → 0.27). Paper's practitioner section and Conclusion
  updated; the referee's "untested regimes" objection is closed.**
  → [experiments/019_claimed-regimes](experiments/019_claimed-regimes/README.md)
- ~~E018 transfer calibration (critique T2.4)~~ **Done (2026-07-04, job
  11677649): plain within-family angle transfer from one n=10 source beats
  every proxy variant on all 28 cells, usually by ~10× (pooled regret 0.007
  vs 0.06; only dense-ER p=3 strains, 0.008 → 0.042 by n=18). Parameter
  concentration made concrete; the paper's abstract and §5.4 now state the
  symmetric practical conclusion — the proxy's regime is where no cheap
  source landscape exists.**
  → [experiments/018_transfer-calibration](experiments/018_transfer-calibration/README.md)
- ~~E017 error directions (July re-test of 012's question)~~ **Done — the
  missing half of the 012 story: with gauge-fixed metrics and controlled error
  directions at matched entrywise MSE, the weighted transfer error does govern
  regret (ρ=0.74; invisible directions are free even at 50% error), shape
  families saturate it (which is why 012 saw no scalar signal), weighted-norm
  refitting does not rescue them, and normalizing the objective is an
  instrument, not a recipe (it hurts the analytical proxy 0.05 → 0.16).**
  See the joint What-we-know entry.
  → [experiments/017_error-directions](experiments/017_error-directions/README.md)
- ~~Open framing question: at p=3 the proxy-chosen schedules equalize leakage, so
  is regret there governed by landscape/argmax robustness rather than fidelity?~~
  **Resolved (exp 014): argmax transfer.** At p=3 fidelity decouples from regret
  (ρ=−0.02) while argmax displacement still predicts it (ρ=0.63); landscape
  flat-peak robustness does not (ρ≈0.1). §5 should narrate depth as a
  parameter-space (argmax-transfer) story, not a fidelity one.
  → [experiments/014_argmax-robustness](experiments/014_argmax-robustness/README.md)
- ~~T2.0 theory write-up~~ **Done** — `theory_compression.tex` (standalone
  proofs) and the paper draft's §3 in the paper repo (see Spencer decision 5).
- Phase-4 assembly progress (2026-07-04, paper commit `e9ba232`, 7 pp,
  compiles clean): number-audit fixes, abstract/intro reframe,
  related-work rewrite with 25 verified references from the union of both
  lit passes, E015 lemmas in §4, E014 argmax-transfer + figure in §5, E016
  in the Discussion, and the exps-012+017 fitted-shape paragraph in §5.4 —
  all done. §5.3 ranking figure + regret table restored (paper commit
  `f10027b` — the generator's outputs were already committed in exp 004; no
  rerun needed). Remaining: optional E015/E016 figures,
  dense-sentence pass (critique T3.1), glossary table, author list,
  acknowledgments, and Spencer decisions 1–6 below.

## Decisions (resolved autonomously 2026-07-04, per Spencer's directive)

Spencer delegated everything except administrativa ("use your best judgment,
don't wait"). Resolutions, auditable in the journal:

1. **Author list + acknowledgments** — parked (Spencer said not to worry
   about these; placeholders remain in the tex).
2. **Venue/format** — the current draft is already
   `\documentclass[conference]{IEEEtran}`; staying with QCE 2027 conference
   format, 8 pp of a ~8–12 pp budget. *(Decided.)*
3. **Headline framing** — mechanism + argmax-transfer lead; the family
   ranking is corroboration with its scope shown (Fig. ranking panel b).
   *(Decided; implemented.)*
4. **E014/E015/E016 fold-in** — done (all three are in the paper with
   figures/lemmas). *(Decided; implemented.)*
5. **`theory_compression.tex`** — retired: header now marks it superseded
   by the paper's §III–IV; kept for proof-detail history, not maintained.
   *(Decided; implemented.)*
6. **§6 framing after exps 012+017** — adopted ("no *absolute* information,
   but the raw landscape is argmax-informative"; direction > size).
   `explainer.md` is the living results walkthrough; `paper_explainer.md`
   kept as a frozen paper companion with a pointer header. *(Decided;
   implemented.)*
