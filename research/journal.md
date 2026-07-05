# Autonomous work journal

*Running log of what the autonomous loop did, decisions taken without Spencer,
and questions queued for Monday. Newest entry first. Established results go in
[STATUS.md](STATUS.md); plain-language versions in [explainer.md](explainer.md);
this file is the working diary.*

## 2026-07-05 (Monday, morning — verification pass on the E025–E031 additions)

### Done

- **Number-verification pass on everything added this morning: 8 findings,
  all applied** (paper `1f7d9eb`; report `verify_pass_2026-07-05b.md`).
  Worst was a genuine self-contradiction ("winning cleanly for α>2" vs
  E026's α=2.2 proxy win eight lines later — now "cleanly at α≥3"); plus
  an honest tie stated as a win (lognormal-depth best source 0.0729 vs
  proxy 0.0731), per-stage evaluation budgets, two roundings, the
  abstract/Conclusion moved to the dispersion framing, and Limits now
  inventories E026 and the 5-instance refinement cells. Everything else —
  including the V₂ Lemma/Prop against the derivation — verified correct
  from the raw CSVs. Five review passes total have now audited the draft.

## 2026-07-05 (Monday, morning — E031: the guarded recipe, validated everywhere)

### Done

- **E031 complete** (array 11782218): the composed recipe with the
  two-evaluation guard — map-init, surrogate ramp-polish, keep the
  better — never hurts (guard recovers the unweighted case where polish
  degrades a ceiling-level init) and gains +0.05 (proxy niche) to
  +0.09–0.12 (hard corner, best-of-two 0.92–0.95). Paper's closing
  recipe updated to the guarded form with per-regime numbers; range
  E001–E031. The morning arc: E028 (last map cell) → V₂ derivation
  (open problem closed) → E029 (corner = grid artifact) → E030 (free
  surrogate polish) → E031 (guarded recipe validated). The program is
  now complete *and* ends on a fully quantified, deployable procedure.

## 2026-07-05 (Monday, morning — E030: the proxy's second job is free polishing)

### Done

- **E030 complete** (array 11779172): compass refinement on the binned
  proxy's own landscape — free — lifts true AR +0.05–0.12 in the hard
  corner within the ramp family (one confirmation eval), while full-2p
  surrogate refinement overfits model error under extreme tails (Pareto:
  0.83 < 0.88). The practitioner section now ends constructively:
  initialize from the map, polish on the surrogate within the schedule
  family, confirm once. Paper + docs updated; range E001–E030.

## 2026-07-05 (Monday, morning — E029: the hard corner was the grid's fault)

### Done

- **E029 complete** (array 11773952): compass refinement from any of four
  grid starts converges to essentially one optimum in both E028 worst
  cells (full-2p residual 0.000–0.006); continuous ramp endpoints alone
  recover most of it; the 8⁴ grid ceilings themselves were 0.02–0.03 low.
  E028's corner reframed in the paper: sharp landscape, not traps —
  local refinement erases it wherever true evaluations are affordable;
  it remains hard only in the zero-refinement (shots-limited) setting.

## 2026-07-05 (Monday, morning — the open math problem falls; E029 running)

### Done

- **The V₂ conditioning correction is derived** (theory agent, machine-
  verified): exact quartic split T = S² + m − R₄ (S = m − 2c) makes
  E[T|c] exactly quadratic plus a small remainder; the quadratic
  projection Vq = (36τ²D − 72τ²W + mW²)/(mD − 36τ²) with W = 4P₁+16q,
  D = 2m(m−1)+24q captures 94–99% of the measured correction on every
  instance. All three mystery constants resolved: 36 = linear regression
  only; measured ≈60 = finite-n Vq; triangle-free residual ≈32 = W²/D;
  ER asymptote κ(p) with κ(½) = 45. In the paper as Lemma 6 + Prop 7
  (`19e8a96`); derivation + verification script committed. The
  Discussion's "open and looks tractable" is now "nearly closed" with
  only the degree-≥3 residue (4–6%) remaining.
- **E029 running** (array 11773952): grid-artifact-vs-intrinsic test of
  E028's hard corner via compass-search refinement.

## 2026-07-05 (Monday, morning — loop resumed on request; E028 closes the map's last cell)

### Done

- Spencer said "continue the research loop" — resumed full pace.
- **E028 complete** (array 11770562): (dispersion × depth), the last
  unmeasured cell. Under lognormal σ=2.5 at p=10/20 the proxy remains the
  most reliable single-shot method (0.05–0.11 vs single sources 0.07–0.29,
  a lottery; pooled sources recover on sparse); under Pareto 1.5 at depth
  every method degrades toward parity. The honest completion: the map now
  also names the corner where nothing tested works well — a marked target
  for future methods, stated as such in the paper.

## 2026-07-05 (Monday, small hours — E027: the map is not a MaxCut story; night closed)

### Done

- **E027 complete** (array 11757023): Max-3-XOR replicates every pillar —
  exactN p=1 regret 0.031–0.040, sampledN ≥ exactN (4th occurrence,
  dramatically so at 8n/16/p=3: 0.044 vs 0.086), transfer/universal
  near-perfect (≈0.000) on unweighted instances, λ ∝ clause density; the
  MaxCut cubic cancellation is confirmed as the only problem-specific
  piece (large-angle λ 0.84–0.95). Paper Limits updated; explainer 26.
- **The loop now rests at an hourly heartbeat.** Every frontier item named
  by any referee pass or by the Limits section has been measured
  (E019–E027, all in one night); what remains is Spencer-facing: the
  length-cut menu, author list, and reading the E022–E027 arc. STATUS and
  this journal are the Monday entry points.

## 2026-07-05 (Monday, small hours — E026: the called shot lands)

### Done

- **E026 complete** (array 11756459): pre-registered prediction (flip in
  α ∈ (2.2, 2.6) on Pareto, from CV² ≈ 1) **confirmed on dense ER** —
  proxy wins at 2.2, exact tie at 2.4, transfer at 2.6. 3-regular keeps
  transfer throughout (its flip is below α=2): the order-unity dispersion
  law holds with a family-dependent constant. Paper sentence + README +
  explainer 25 + STATUS updated. The E022–E026 arc is now: niche found →
  boundary located → mechanism identified → prediction confirmed.

## 2026-07-05 (Monday, small hours — E025: the mechanism is dispersion)

### Done

- **E025 complete** (array 11755955, ~1.5 min/task): lognormal weights
  reproduce the crossover with variance finite at every σ — flip at
  σ≈1.0–1.5 (CV² ≈ 2–8), transfer collapsing to 0.16–0.35 by σ=2.5 while
  the proxy improves to 0.004. E023's "infinite-variance boundary" was a
  grid artifact (α-grid jumps CV² from 0.33 to ∞ between α=3 and 2). The
  paper's quotable claim upgraded and made more practical: **the
  transfer↔proxy boundary is order-unity weight dispersion (CV²), which a
  practitioner can compute from the weights in one line**; infinite
  variance is sufficient, not necessary. Paper `8337edd`; explainer 24;
  Limits updated (lognormal no longer "untested").
- Follow-up left open (journaled, not run): a targeted Pareto α∈(2,3)
  mini-sweep to confirm the CV²-threshold unification on the Pareto side
  (predicted flip near α≈2.2–2.4).

## 2026-07-05 (Sunday, night, last — final referee verdict applied: weak accept

### Done

- **Final referee pass returned: weak accept, mechanical fixes only** —
  and all six defects are applied (paper `a73413b`, now 10 pp): the
  falsified Limits paragraph rewritten with the honest weighted scope +
  the lognormal and Max-k-XOR mitigations; "locates the boundary
  precisely" softened to "consistent with α≈2" (the reviewer recomputed
  the CSVs: 3-regular already flips at α=2 at >4σ — sharp point is
  dense-ER-only); instance counts and SE scoping added to every E022–E024
  claim (two headline cells are 2×SE ties and now say so); the proxy's
  α-range quoted honestly (0.012–0.046); abstract's 2–7× gains its p=1
  qualifier; contribution 5 gains the boundary. Review on file:
  `research/final_referee_2026-07-05.md`.
- **Not applied (deliberate, queued for Spencer):** the reviewer's
  length-cut menu (drop fig:bound's right panel, fold fig:ranking to one
  panel, compress the V₂ lemmas, single-column fig:binned) — the draft
  fits 10 pp as-is, and those cuts trade real content for margin; a call
  for Monday.

### Night summary (for Monday's read)

Since Spencer's "never idle" directive: experiments E019–E024 designed,
run, and folded in (claimed-regimes test; weighted binning; weighted
transfer rematch; heavy-tail crossover; α-sweep boundary; scale
confirmation), three referee/coherence passes applied, and the paper went
from an honest-negative story to a complete practitioner map with a
measured positive niche at a principled boundary. Everything committed,
pushed, reproducible; verdict on file: weak accept pending human polish.

## 2026-07-05 (Sunday, night — E024: boundary scale-confirmed; final referee pass launched)

### Done

- **E024 complete** (GPU array 11754145): at n=22/26 with the sampled-only
  binned proxy, α=1.5 keeps the proxy ahead in every cell (0.015–0.025 vs
  0.030–0.130 all transfer sources) and α=3 keeps transfer ahead in every
  cell. The infinite-variance rule is scale-robust and delivered by the
  practical method. Paper sentence added; explainer result 23; STATUS.
- **Final full referee pass launched** on the completed 9-pp draft
  (E022–E024 additions included) — last substantive quality gate before
  Spencer's Monday read.

## 2026-07-05 (Sunday, night — E023: the boundary is the infinite-variance point)

- **α-crossover figure added** (`abc0127`): proxy flat, transfer sweeping
  through it, crossing on the α=2 line — the practitioner section's
  punchline visual. Next queued: E024, confirming the α≤2 proxy win
  survives scale (n=22/26, GPU ceilings, Pareto 1.5 vs 3.0).

### Done

- **E023 complete** (array 11751921, 10 tasks, ~2.5 min each): sweeping
  Pareto α ∈ {1.2, 1.5, 2, 3, 5} shows the proxy's regret nearly
  α-independent (~0.02–0.03 at p=1) while transfer moves through it —
  cleanly winning for α>2, tying at α=2 (dense ER: 0.028 vs 0.028),
  losing below. **The crossover coincides with the weight variance
  becoming infinite** — a principled boundary that pre-empts the
  "cherry-picked α" objection. At p=3 the pattern shifts one rung heavier
  and the transfer source-lottery becomes extreme (0.035–0.212 from
  same-family sources at α=1.2). Paper updated (`5302abc`); explainer
  result 22; STATUS.

## 2026-07-05 (Sunday, later — E022: the crossover found; practitioner map complete)

### Done

- **E022 complete in minutes** (array 11751040, 8 tasks). **The headline
  the whole calibration arc was missing: under Pareto(α=1.5) weights the
  binned per-instance proxy WINS** — every p=1 cell by 2–7×; stable at
  p=3 (0.046–0.067) against lottery-like transfer sources (0.021–0.140).
  Exp(1) remains transfer territory, confirming the heterogeneity dial:
  U[0,1] < Exp(1) ≪ Pareto(1.5). Mechanism: mean-weight rescaling fails
  when the mean is tail-dominated; quantile bins absorb giant edges per
  instance. Third independent observation of sampled-N ≥ exact-N.
- **Paper updated end-to-end** (`7613c37`, 9 pp clean): abstract,
  contribution 5, the practitioner section (E021 paragraph now flows into
  the E022 inversion), and a Conclusion that names the completed map —
  family concentration ⇒ transfer; rescaling-defeating heterogeneity ⇒
  per-instance binned proxy; the calculus tells you which regime you're
  in. Explainer result 21; STATUS updated.
- With E022, the story arc is genuinely complete: theory (Thms 1–3 +
  V₂ anatomy + binning extension), anatomy (density law, argmax
  transfer), honest calibrations (E013/018/019/021), and a real positive
  niche (E020+E022). Twenty-two experiments, all reproducible.

### Open questions for Monday (standing)

- Author list + acknowledgments (parked per Spencer).
- Optional next experiments if the loop continues: E015/E016 optional
  figures; a Pareto-α sweep to locate the crossover boundary; n=22–26
  heavy-tail confirmation (GPU); a fresh full referee pass on the final
  9-pp draft.

## 2026-07-05 (Sunday — coherence pass applied; E022 queued)

### Done

- **Coherence pass: 16 findings, all resolved** (paper `d2c0c6c`; report
  committed). Worst: Table II's caption still carried the retracted
  "regret is flat in n" claim directly above a table showing monotone
  growth — exactly the kind of survivor a hostile referee screenshots.
  Also fixed: three wrong-target section refs after the practitioner
  split; "best proxy variant" scoped to p≤3 everywhere; hardware shots
  removed from the timing regimes (transfer avoids them equally — the
  argument only separates offline from on-device methods); η_F jargon;
  SE conventions split (AR means 0.002–0.006 vs regret cells ~0.001);
  bin-count phrasing aligned; Table III's K renamed J (with the n factor
  restored in the sweep scaling); the weighted extension added to
  contribution 5.

### Next (E022, the finish line for the practitioner story)

- **Strong heterogeneity with a FAIR transfer baseline.** E021's i.i.d.
  U[0,1] weights weren't heterogeneous enough. E022: heavy-tailed
  per-edge weights (Pareto α=1.5; Exp(1) as the milder rung) where the
  landscape shape itself varies between instances — and, critically, give
  transfer its practitioner's fix (γ rescaled by the source/target mean
  weight ratio) so the proxy only wins if per-instance structure beyond
  scale matters. Same cells/machinery as E020/E021. If rescaled transfer
  still wins, the practitioner section closes with "always transfer";
  if the binned proxy wins, we finally have its niche, measured.

## 2026-07-05 (Sunday, early — E021: the weighted rematch; coherence pass running)

### Done

- **E021 complete in minutes** (array 11748691; reused E020 targets and
  the allocation-free kernel). Weights erode transfer 3–10× and destroy
  the universal schedule on dense ER at depth, but same-family weighted
  sources still win every cell against binned exact/sampled proxies.
  Direction (heterogeneity → per-instance methods) confirmed; crossover
  not reached at i.i.d. U[0,1]. Folded into the paper (`9beec01`),
  explainer result 20, STATUS.
- **Next open question made precise:** the proxy's viable niche, if any,
  is strong weight heterogeneity (heavy tails, mixed scales, structured
  weights) — a concrete follow-up experiment, not a vague hope.
- Coherence-only referee pass on the heavily-edited draft is running in
  the background; findings get applied when it reports.

## 2026-07-05 (Sunday, early — E020 lands: weighted MaxCut via binning works)

### Done

- **E020 complete** (tasks 1/2/4 in job 11744099; task 3 rerun as 11745830
  after two Julia-1.12 GC aborts — root-caused to allocation churn in the
  threaded sweep, fixed by memoized factors + per-thread buffers,
  bit-identical smoke, and 3× faster). **Headline: the framework's
  integer-cost "wall" is gone.** λ² ≈ λ_struct² + O(1/K) measured across
  both families and sizes; binned-proxy regret saturates by K≈64–128 at
  the unweighted exact-compression levels (ER(0.5) n=16 p=1: 0.038 vs
  ~0.032; 3-regular n=16 p=3: 0.103 vs 0.105). Folded into the paper
  (`1d90d9b`): practitioner-section paragraph, abstract clause, Limits
  rescoped from "not analyzed" to measured scope. Explainer result 19.
- **Ops notes:** two core dumps from the GC crashes were nearly committed
  by a background `git add` (multi-GB, slow hash over NFS — looked like a
  hang; an index.lock removal race followed). Cleaned up: core dumps
  deleted, `core.*` gitignored, waiters now use sacct state checks instead
  of squeue emptiness.

## 2026-07-04 (Saturday, night, latest — E019 lands: both claimed regimes go to transfer)

### Done

- **E019 complete and folded in** (array 11740470 all COMPLETED; paper
  `c346ad2`; STATUS + explainer result 18 updated). Depth p=10/20: every
  proxy variant collapses (exact 0.07–0.14, sampled 0.11–0.21, analytical
  0.18–0.37) while transfer reaches 0.001–0.03 — the exact compression
  misplaces the argmax by ~0.1 AR at fidelities of 0.7–0.8, the cleanest
  argmax-transfer demonstration we have. Beyond the exact-N wall (n=22/26):
  transfer wins every cell; dense-ER p=3 worsens with n (0.047→0.061) but
  stays ahead of sampled (0.078–0.081) and analytical (0.083–0.106); the
  pooled universal schedule beats the proxies almost everywhere. Positives
  worth keeping: the analytical argmax improves with n on sparse regular
  graphs (asymptotic derivation, now measured) and its binomial P is
  load-bearing (empirical P degrades 3-regular p=3 from 0.05 to 0.27).
- The referee's #2 risk ("advertised regimes untested") is closed with
  data; the p=3 depth gap is closed (comparisons now span p=1–20,
  n=12–26).

### Next experiment (E020, designing now)

- **Binned compression for weighted MaxCut** — converts the referee's #1
  risk (integer costs are load-bearing; singleton classes trivialize the
  compression) into a new contribution. Random U[0,1] edge weights make
  every cost class a singleton; quantile-bin the costs into K integer
  labels and drive the existing machinery with them (projector = bin
  averaging; proxy phases use bin-mean costs; true evolution uses true
  weights). Measure per-layer leakage vs K (predicted: compression term +
  binning term shrinking in K) and binned-proxy parameter-setting regret
  vs true ceilings, n=14–16, p=1 and p=3. Theorems 2–3 apply to any fixed
  partition, so the calculus extends verbatim — the question is whether
  the binning term is small enough at modest K for the proxy to remain
  usable, and how the K-tradeoff looks.

## 2026-07-04 (Saturday, night, later — fresh referee pass applied)

### Done

- **Fresh referee pass on the current draft came back accept-trajectory
  with one real factual error**, now fixed (paper `f3f62f0`): the bolded
  "regret does not grow with n" contradicted the paper's own per-family
  regret table (monotone growth in all 7 families, +0.005 to +0.018 over
  n=12→18, resolvable against the ~0.001 SEs). The text now states the
  numbers and owns the correction. Report:
  `research/referee_pass_2026-07-04.md`.
- **Overselling residues cleaned**: sampled-N demoted from "recommended
  recipe" to "best proxy variant" with an explicit pointer that transfer
  owns the family regime; the hardware claim rescoped (transfer's source
  is classical too — shots don't rescue the proxy; the honest
  differentiator is per-instance vs family-dependent methods); "within
  0.01" → measured 0.011–0.016; "decouples entirely" scoped to p=3;
  GPU-parity scoped to n≤20 with the extrapolated n≈22–24 crossover.
- **Limits expanded** with the weighted-cost boundary (singleton classes
  trivialize the compression; binned classes flagged as the open
  extension; Thms 2–3 are partition-agnostic, the cubic lemma is
  MaxCut-specific) and the noiseless scope of all comparisons.
- **Structure**: abstract cut ~330 → ~200 words; "What should a
  practitioner do?" split into its own subsection; a short Conclusion
  added; the stale "literature pass predates submission" sentence
  removed. Skipped by judgment: compressing the V₂ lemmas and merging the
  two trajectory passages (both would cost clarity more than space).
- Of the referee's top-5 unanswered questions, three are now answered in
  Limits (weighted costs, noise, MaxCut-specificity), the empty-toolbox
  point is in the practitioner subsection + Conclusion, and the p=3 depth
  gap is exactly what E019 (running, array 11740470) closes.

## 2026-07-04 (Saturday, night — Spencer's new directive: never idle, decide and go)

### Directive (from Spencer, verbatim intent)

Skip administrativa (author list etc.); decide judgment calls (wording,
framing) with best judgment; each loop iteration must critically review the
paper (methodology, novelty/publishability), plan and run new experiments,
and keep improving — never stop to wait for input. Saved to persistent
memory so future sessions start in this mode.

### Done immediately

- **E019 designed, smoke-tested, submitted (GPU array 11740470, 10 tasks).**
  The paper's biggest exposure after E013+E018: the claimed-value regimes
  were untested. Part A pushes depth (p=10, 20 ramps; ER(0.5)/ER(0.25)/
  3-regular; n=14–16) where transfer showed its only strain; Part B goes
  beyond the exact-N wall (n=22, 26 — exact N impossible, ceilings still
  GPU-computable) where only analytical/sampled N exist. Methods: ceiling,
  exactN (A only), sampledN, paper (binP + empP), transfer×3, and a pooled
  "universal" schedule. Either the proxy earns a regime or the honest
  conclusion hardens — both are publishable outcomes.
- **Fresh adversarial referee pass launched** (background agent) on the
  CURRENT 8-pp draft — the June critique reviewed a longer, older version.
  Focus: methodology soundness, novelty vs the newly cited neighbors,
  coherence after the E017/E018 additions (any sentence still overselling
  the proxy), structure, and the top unanswered referee questions. Report →
  `research/referee_pass_2026-07-04.md`.
- **Decisions 2–6 resolved autonomously** (see STATUS): QCE/IEEEtran stays;
  mechanism-led framing stays; snippets already folded;
  `theory_compression.tex` retired with a superseded header;
  §6 refined claim adopted; `explainer.md` designated the living doc with a
  pointer header added to `paper_explainer.md`. Author list parked per
  Spencer.

## 2026-07-04 (Saturday, evening — E018 run and folded in; pre-Spencer backlog complete)

### Done

- **E018 complete in 65 s** (job 11677649; 700 targets regenerated from the
  002/010 seeds, all asserted): **transfer wins every cell** — pooled
  regret 0.007 (median 0.003) vs 0.06 for the exact-compression proxy;
  best-of-3 sources drives 3-regular to ~0.0003 and WS(0.1) to ~0.0001;
  the only strain is dense ER at p=3 (0.008 → 0.042 as n grows to 18),
  the same regime where compression leaks most.
  → `experiments/018_transfer-calibration/` (+ analysis.txt)
- **Folded into the paper immediately** (`e7dea19`, now 8 pp): a
  calibration paragraph in §5.4 citing the concentration/fixed-angle
  literature, and one honest closing sentence in the abstract ("the
  proxy's regime is where no cheap source landscape exists"). Explainer
  gained result 17; STATUS updated.
- **Noticed while editing:** the current draft header is already
  `\documentclass[conference]{IEEEtran}` — the critique's venue-ambiguity
  item (decision 2) referred to the lost long draft's header; what's left
  of that decision is only the page-budget question (8 pp now, QCE limit
  ~8–12).

### State of the loop

With E018 in, **every item actionable without Spencer is done**: critique
burn-down complete, both lit passes folded, all experiment results written
into paper + explainer + STATUS, both repos pushed and in sync. What
remains is Spencer-gated (decisions 1–6, §6/§5.4 wording sign-off, page
budget) or optional (E015/E016 figures pending venue, a T2.2 ruggedness
mechanism experiment, deeper readability sweeps). The loop will keep
ticking at a slower cadence and pick up anything new (e.g., replies,
CI, or a fresh idea vetted against the program plan).

## 2026-07-04 (Saturday, later — critique burn-down complete; E018 queued)

### Done

- **§5.3 ranking figure + regret table restored** (paper `f10027b`) — the
  last unrecovered June-13 item; exp 004 had the outputs committed all
  along, so this was pure wiring. Captions carry the honest scope (panel b
  shows the p=3 collapse) and distinguish proxy-chosen-angle leakage from
  the fixed-angle density law.
- **Abstract's four-clause closer split; argmax transfer named up front**
  (`ad09027`); **model-error subsection opens with its one-line thesis;
  glossary reading key rebuilt** and referenced from the intro
  (`ad50808`). With these, every critique item that doesn't require
  Spencer is closed: T1.1, T1.2, T1.4, T2.1, T2.3, T2.5, T3.1 (named
  items), T3.3, T3.5 done; T1.3/T3.2 are decisions 2/5; T3.4 has no
  target section in this draft.
- **Explainer: plain-language sections 14–16** for the June-line
  experiments (argmax transfer at depth; the V₂ codegree/triangle
  anatomy; the cheap-prefix-frame obstacle) (`a24abf7`).

### Next

- **E018 (critique T2.4): the transfer-calibration row.** Question: how
  does proxy regret compare to plain angle transfer from a small instance
  of the same family? Design sketch: per (family, n) cell of exps
  002/010, grid-optimize angles on one n=10 source instance's true
  landscape, apply to every target instance, and measure regret against
  the targets' already-recorded grid ceilings (read from the 002/010
  CSVs; seeds are in the CSVs, so instances regenerate exactly). One
  statevector evaluation per target per depth — light, but routed through
  Slurm per the ground rules. Produces the one-row calibration a referee
  will ask for ("transfer regret X vs proxy regret Y on the same cells").

## 2026-07-04 (Saturday — paper session: recovered edits + both lit passes folded in)

### Done

- **The paper absorbed everything the reconciliation staged** (paper repo
  commit `e9ba232`, pushed; compiles clean at 7 pp, zero unresolved refs):
  - Number-audit fixes re-applied (the −0.007 misattribution, the
    overlap-range phrasing); experiment range bumped to E001–E017.
  - Abstract/intro reframed per critique T2.3: no more
    "depths where simulation is intractable" promise; an explicit
    "not a laptop speedup" expectation-setting paragraph in the intro.
  - Related work rewritten from the *union* of the June and July lit
    passes: 25 verified references added (Grover-mixer = exact
    homogeneity; Kemeny–Snell lumping origins; Burgholzer bisimulation
    rescopes "to our knowledge"; Dirac–Frenkel lineage for Thm 2 with
    Lubich/Martinazzo/Zoufal, claimed as known-technique-new-object;
    the He et al. parameter-setting toolbox; the classical-surrogate
    family; Krüger–Mauerer contrast).
  - E015's two V₂ lemmas + density-law remark now in §4; the
    linear-depth law's proved-vs-measured status stated plainly (T2.5);
    the |+⟩-eigenstate observation completes the cubic sketch (T3.3).
  - E014 argmax-transfer paragraph + figure in §5; E016 resolves the
    Discussion's cheap-frame question.
  - New §5.4 paragraph: the fitted-shape paradox dissected via exps
    012+017, including the refined "no *absolute* information, but the
    raw landscape is argmax-informative" claim and the
    normalization-as-instrument finding.
- **Remaining paper backlog:** §5.3 ranking figure + regret table (the one
  still-unrecovered June-13 item; generator committed in exp 004 — needs a
  Julia/CairoMakie run, queued next); optional v2_anatomy and prefix_frame
  figures; T3.1 dense-sentence pass; T3.5 glossary table (existed only in
  the lost long draft); Spencer decisions 1–6 (venue, authors, framing,
  theory_compression.tex, §6 wording).

## 2026-07-03 (Friday, evening — branch reconciliation + second lit scan)

### The discovery

Pushing today's work failed: `origin/ClaudeResearch` held **21 commits this
clone had never fetched** — a June-13 overnight autonomous session
(timestamps +0900) that produced experiments **012 (fitted-shape paradox),
014 (argmax-robustness), 015 (V₂ anatomy), 016 (cheap-prefix frame)**, the
systematic literature pass, the paper critique with number audit, and
ready-to-paste paper snippets. The two lines forked at `8680853`. This
clone's July sessions unknowingly duplicated parts of it (the E013 write-up
describes the *same* Slurm runs; the kernel fix trees are byte-identical)
and re-tested exp 012's question under the name "E014."

### Done

- **Merged `origin/ClaudeResearch`** (merge commit `9c49d7b`): STATUS.md
  rewritten to weave both accounts with a prominent reconciliation note; the
  two E013 READMEs (same run, two write-ups) unified.
- **Renamed today's experiment 014 → `017_error-directions`** (the June line
  owns 012 and 014; pushed history wins naming). Slurm job name/log
  identifiers keep `e014`/11575420 in sacct; a README note maps them.
- **Wrote the 012 ↔ 017 synthesis** (new README section + joint STATUS
  entry): 012 found *no scalar norm* discriminates regret across its
  fitted-model zoo (only argmax displacement, ρ≈0.7); 017 shows the weighted
  transfer error *is* causal when given dynamic range (matched-MSE direction
  quartet) and that shape fits *saturate* it — which is exactly why 012 saw
  no scalar signal. The two designs interlock rather than contradict.
- **Determined the June-13 paper-repo edits are lost** (committed only on an
  unpushed clone: §5.3 figure inclusion, three draft fixes, two number
  corrections, five citations; the critique references a ~800-line tex, ours
  is ~550). All content is reconstructible from committed sources
  (`paper_critique.md`, `literature_pass.md`, figure generators,
  `proposed_paper_additions.md`) — queued as the next paper-thread task.
- **Second literature scan complete** (`research/lit_scan_2026-07-03.md`) —
  run independently before the June pass was discovered, so it doubles as a
  cross-check. Agreements: novelty of the four core contributions confirmed
  by two independent scans; Krüger–Mauerer (Quantum 2025) is the closest
  neighbor in both. **New referee risks the June pass missed:** (i)
  Theorem 2's technique is the *Dirac–Frenkel a posteriori bound* (Lubich
  2008, Thm II.1.5), already ported to variational quantum evolution by
  Zoufal et al. 2023 — cite and reframe as "known technique, new object";
  (ii) Burgholzer et al. (ACM TQC 2025) apply lumpability-style bisimulation
  to quantum circuits including QAOA-on-MaxCut (exact-only) — the "to our
  knowledge" sentence needs rescoping; (iii) Grover-mixer QAOA is *exactly*
  Perfect Homogeneity — a zero-leakage boundary case Theorem 1 should
  mention; (iv) the draft's parameter-setting paragraph needs ~7 standard
  refs (Zhou FOURIER/INTERP, TQA, fixed angles, concentration, linear-ramp
  at scale); (v) citation hygiene: the tex doesn't load
  `Mitsubishi_B_references.bib` at all (inline thebibliography), sud2024
  title is the arXiv variant, a 2026 "subspace compression" sketching
  preprint warrants a disambiguating footnote.

### Decisions taken autonomously

- **Merged rather than rebased** (preserves both true histories), resolved
  doc conflicts by weaving, renumbered my experiment rather than theirs
  (theirs is pushed), and recorded everything here + STATUS for Monday.
- **Both lit-pass reports stay**; §2 revision will use their union
  (June: breadth + citer sweep; July: technique-lineage risks).
- `explainer.md` (results walkthrough) and `paper_explainer.md` (paper
  companion) both kept for now — consolidation queued as Spencer decision 6.

### Open questions for Monday (added)

8. **Why did this clone never fetch June 13–July 3?** Worth a workflow
   guard: the loop now runs `git fetch && git status -sb` at session start
   (added to the loop's own checklist) so divergence surfaces immediately,
   not at push time.
9. **The June-13 session's paper edits are gone** — re-applying them is
   mechanical from the critique + generators, but if you have that clone
   (another machine/directory?), a `git push` from it would save an hour and
   preserve exact history. Check before Monday's session if convenient.

## 2026-07-03 (Friday, afternoon — continuous-loop session)

### Done

- **E017 COMPLETE, same afternoon** — all 10 array tasks finished in 2–8 min
  each (far under the 90-min budget), 3450 rows, logs clean, analysis in
  `experiments/017_error-directions/analysis.txt`. Results (details in
  the experiment README):
  1. *Paradox resolved:* weighted transfer error governs regret (pooled
     Spearman 0.74 vs 0.38 for MSE); invisible perturbations at 50%
     entrywise error carry zero excess regret; MSE orders triangle-vs-normal
     *against* regret in 105/150 instances (an anti-predictor).
  2. *Prescriptive fix fails:* weighted-norm refits don't rescue the shapes
     (triangle +0.007 marginal, normal −0.016; representational error stays
     ≥ 0.92). Shape fitting is dead on representational, not
     objective-function, grounds.
  3. *Normalized objective is instrument, not recipe:* it hurts exact N
     slightly (0.0315 → 0.0415) and the analytical proxy badly
     (0.05 → 0.16, worse in ~132/150) — PaperProxy's raw landscape is
     argmax-informative, refining the 004–006 "values are noise" claim.
- **E017 implemented and submitted: Slurm array job 11575420** (10 tasks =
  5 families × n ∈ {12,14}, 15 instances each, CPU partition since n ≤ 14
  needs no GPU). Three parts as designed this morning: (A) matched-MSE
  perturbations with four d-profiles, (B) triangle/normal entrywise-MSE fits
  + PaperProxy, (C) the same fits in the Theorem-3 weighted norm.
  → `experiments/017_error-directions/`
- **The smoke run caught a real methodological discovery, not just a bug:**
  under the paper's raw Eq.-9 objective (exps 001–010 convention), *any*
  coherent perturbation of N — even ε = 0.01 and transfer-invisible — plants
  a norm-inflated beacon (predicted ⟨C⟩ ~ ε², up to 5×10⁵ on a 21-edge
  graph) that hijacks the argmax. This is the exp-004–006 pathology in its
  general form: the raw objective is unusable for model-error studies.
  E017 therefore records a second, **normalized** ranking (divide by the
  compressed state's weight Σ2ⁿP|Q|²; gauge-invariant) for every variant.
- **Smoke-scale confirmation of the Theorem-3 prediction** (one n=10
  ER(0.5) instance): normalized regret at matched entrywise MSE spans
  0.0000 (nulled profile, invisible by construction) to 0.32 (aligned
  profile) — the weighted error, not the MSE, governs regret. Also observed:
  exact-N normw = 0.978 ≤ 1 (Thm-1 contractivity visible in the wild);
  normalization even improved exact-N's own parameter choice; but it *hurt*
  PaperProxy on that instance (raw regret 0.0011 → nrm 0.0995) — watch this
  at scale, it complicates a clean "normalize, don't veto" recipe.
- **Paper repo pushed** (pipeline-cost subsection with E013 timing table,
  ff0e743, was sitting unpushed).
- **Systematic literature pass launched** (background agent: parameter
  setting/transfer, mean-field AOA & classical surrogates, Zwanzig–Mori /
  lumpability precedents, non-ER instance dependence, DOS-based analyses,
  and direct follow-ups to Sud et al.). Report will land at
  `research/lit_scan_2026-07-03.md`; findings get folded into §2 next.

### Decisions taken autonomously

- **E017 ranks parameters under BOTH the raw Eq.-9 objective and the
  normalized one**, reporting both regrets. Rationale: raw keeps continuity
  with exps 001–010; normalized is the only instrument with dynamic range
  for Part A (and doubles as the "normalize, don't veto" test). If the
  normalized recipe survives the full run, it belongs in §6 and possibly
  changes the paper's recommended objective — flag for Spencer.
- **E017 runs CPU-only** (16 threads/task, general-short): at n ≤ 14 the
  statevector grids and homodists are cheap; the GPU queue is the scarce
  resource and E013 already showed the GPU adds nothing at this scale.
- Perturbation amplitudes are relative (∝ ‖N[c',:,c]‖ per slice) and may
  make N entries negative — accepted deliberately; the proxy iteration is
  linear and Part A stress-tests the *metric*, not a physical model class.

### Open questions for Monday (added)

5. ~~Adopt the normalized objective paper-wide?~~ **Answered by the full run:
   no.** The normalized objective is the correct *measuring instrument* for
   model-error effects (Parts A–C all use it) but a *worse parameter-setter*
   (hurts exact N in 99/150, hurts PaperProxy in ~132/150). Plan: Eq. 9
   stays; §6 gains the E017 story — paradox resolution, shape-fitting
   post-mortem, and the refined claim "raw values carry no *absolute*
   information but the raw landscape is argmax-informative." Sanity-check
   this framing with Spencer before it goes in the paper.
6. ~~PaperProxy worse under normalization — smoke fluke?~~ **Confirmed at
   scale** (see 5). The mechanistic *why* — what about the analytical N
   makes its norm inflation correlate with parameter quality — is open and
   would make a nice §6 paragraph or follow-up note. Candidate story: the
   analytical model's inflation is largest where the true landscape is also
   large (both are driven by constructive small-d interference), except on
   dense ER where the multinomial tail misfires. Not yet tested.
7. ~~E017's Part-A device as a paper figure~~ **Done (2026-07-04):**
   `fitted_shape_directions.png` generated (exp 017), included in §5.4
   (paper `6777c36`) — invisible directions ride the exact-N baseline to
   ε=0.5, visible ones leave it at ε=0.01.

## 2026-07-03 (Friday)

### Done

- **E013 (timing benchmark) written up and committed.** Both Slurm logs
  analyzed: job 9807239 crashed at n=20 on a CUDA grid-dimension overflow;
  job 9807549 (after the kernel fix) completed n=12–20. Headline: obtaining
  exact N is the pipeline bottleneck (O(4^n)); sampled N is 50× cheaper at
  n=20; at simulable sizes GPU brute force matches the proxy's speed, so the
  paper's practical section must make the asymptotic/hardware-cost argument
  honestly. → `experiments/013_timing-benchmark/README.md`
- **Committed the pending CUDA fix** (grid-stride homodist kernel, validated
  GPU-vs-CPU at n=10,12) and the `.gitignore` entry for the nested paper repo.
- **STATUS.md refreshed** (was dated 2026-06-10 and still listed T2.0/E2.1 as
  "next up" though both are done).
- **Created `research/explainer.md`** — plain-language walkthrough of all
  results 001–013, to be updated with every new result (this was missing;
  STATUS.md is a state summary, not an explanation).
- **Launched a 4-lens review** of `qce2027_paper.tex` (theory correctness,
  claims-vs-evidence, writing/venue fit, related-work coverage) with
  adversarial verification of factual findings; fixes to be applied to the
  paper repo when it reports back.
- **E017 (fitted-shape paradox, was "E2.4") in design** — see below.

### Decisions taken autonomously

- **Experiment numbering: next experiment is 014, leaving 012 as a permanent
  gap.** There is no `012_*` directory and no trace of one in git history;
  reusing 012 after 013 would make directory order disagree with execution
  order. *[Superseded same day by the branch reconciliation: 012 and 014 DID
  exist on origin/ClaudeResearch, unfetched; the experiment described below
  was renumbered 017. See the evening entry.]*
- **E017 design: two-part test of the Theorem-3 explanation for the
  fitted-shape paradox.** Part A, controlled perturbations of the exact
  empirical N at matched entrywise MSE but very different f_d(β)-weighted
  error (prediction: regret tracks the weighted error, not the MSE). Part B,
  actual TriangleProxy/NormalProxy fits (entrywise-MSE objective, as G-RIPS
  did) plus PaperProxy, checking whether the weighted model error ranks the
  models' parameter-setting regret where entrywise MSE misranks them.
  Fitting is done in Julia with simple random search (the Python
  `fit_proxy_to_real` is not batch-friendly); n=12–14, five families,
  p=1 grid + ceilings as in exps 002/010.
- The weighted model-error metric is defined as the Frobenius norm of the
  one-layer transfer-matrix difference: E_w(β)² = Σ_{c',c} |Σ_d f_d(β)
  ΔN(c';d,c)|² (the γ phases drop out — the metric depends only on β).
  This is exactly the operator the proxy update applies, so it is the
  Theorem-3-consistent notion of model error.

### Open questions for Monday

1. **Framing check (§6 / practical):** are you comfortable with the honest
   E013 conclusion — "at simulable sizes the proxy is not faster than GPU
   brute force; its value is asymptotic + expensive-evaluation settings"?
   It contradicts the original speedup pitch but matches program.md's
   assessment and is defensible.
2. **E017 scope:** Part B fits triangle/normal with entrywise-MSE random
   search (mirroring what G-RIPS did). Should we *also* fit in the weighted
   norm to demonstrate the prescriptive fix ("fit in the right norm"), or
   keep that for a follow-up? Currently planning to include it — it turns a
   diagnosis into a recipe, which strengthens §6.
3. **Timing table in the paper:** planned as a small table in the practical
   section (n=20 row + scaling statement). OK?
4. **Paper review findings** will be summarized here once the review
   completes; high-priority fixes get applied directly, judgment calls get
   queued for you.
