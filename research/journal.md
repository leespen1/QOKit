# Autonomous work journal

*Running log of what the autonomous loop did, decisions taken without Spencer,
and questions queued for Monday. Newest entry first. Established results go in
[STATUS.md](STATUS.md); plain-language versions in [explainer.md](explainer.md);
this file is the working diary.*

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
