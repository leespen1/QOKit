# Autonomous work journal

*Running log of what the autonomous loop did, decisions taken without Spencer,
and questions queued for Monday. Newest entry first. Established results go in
[STATUS.md](STATUS.md); plain-language versions in [explainer.md](explainer.md);
this file is the working diary.*

## 2026-07-03 (Friday, afternoon — continuous-loop session)

### Done

- **E014 COMPLETE, same afternoon** — all 10 array tasks finished in 2–8 min
  each (far under the 90-min budget), 3450 rows, logs clean, analysis in
  `experiments/014_fitted-shape-paradox/analysis.txt`. Results (details in
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
- **E014 implemented and submitted: Slurm array job 11575420** (10 tasks =
  5 families × n ∈ {12,14}, 15 instances each, CPU partition since n ≤ 14
  needs no GPU). Three parts as designed this morning: (A) matched-MSE
  perturbations with four d-profiles, (B) triangle/normal entrywise-MSE fits
  + PaperProxy, (C) the same fits in the Theorem-3 weighted norm.
  → `experiments/014_fitted-shape-paradox/`
- **The smoke run caught a real methodological discovery, not just a bug:**
  under the paper's raw Eq.-9 objective (exps 001–010 convention), *any*
  coherent perturbation of N — even ε = 0.01 and transfer-invisible — plants
  a norm-inflated beacon (predicted ⟨C⟩ ~ ε², up to 5×10⁵ on a 21-edge
  graph) that hijacks the argmax. This is the exp-004–006 pathology in its
  general form: the raw objective is unusable for model-error studies.
  E014 therefore records a second, **normalized** ranking (divide by the
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

- **E014 ranks parameters under BOTH the raw Eq.-9 objective and the
  normalized one**, reporting both regrets. Rationale: raw keeps continuity
  with exps 001–010; normalized is the only instrument with dynamic range
  for Part A (and doubles as the "normalize, don't veto" test). If the
  normalized recipe survives the full run, it belongs in §6 and possibly
  changes the paper's recommended objective — flag for Spencer.
- **E014 runs CPU-only** (16 threads/task, general-short): at n ≤ 14 the
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
   stays; §6 gains the E014 story — paradox resolution, shape-fitting
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
7. **E014's Part-A device (matched-MSE perturbation quartet) could become a
   paper figure** — regret vs ε for the four profiles is the visual proof
   that entrywise size doesn't matter and direction does. Draft next
   session; needs no new compute.

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
- **E014 (fitted-shape paradox, was "E2.4") in design** — see below.

### Decisions taken autonomously

- **Experiment numbering: next experiment is 014, leaving 012 as a permanent
  gap.** There is no `012_*` directory and no trace of one in git history;
  reusing 012 after 013 would make directory order disagree with execution
  order.
- **E014 design: two-part test of the Theorem-3 explanation for the
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
2. **E014 scope:** Part B fits triangle/normal with entrywise-MSE random
   search (mirroring what G-RIPS did). Should we *also* fit in the weighted
   norm to demonstrate the prescriptive fix ("fit in the right norm"), or
   keep that for a follow-up? Currently planning to include it — it turns a
   diagnosis into a recipe, which strengthens §6.
3. **Timing table in the paper:** planned as a small table in the practical
   section (n=20 row + scaling statement). OK?
4. **Paper review findings** will be summarized here once the review
   completes; high-priority fixes get applied directly, judgment calls get
   queued for you.
