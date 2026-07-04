# Autonomous work journal

*Running log of what the autonomous loop did, decisions taken without Spencer,
and questions queued for Monday. Newest entry first. Established results go in
[STATUS.md](STATUS.md); plain-language versions in [explainer.md](explainer.md);
this file is the working diary.*

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
