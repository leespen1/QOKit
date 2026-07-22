# 24-hour improvement loop journal

*Loop started 2026-07-22 18:26 JST; hard stop after 2026-07-23 18:26 JST.
Spencer's directives: (1) hunt and fix paper errors; (2) make
`paper_explainer.md` complete enough to rewrite the paper from; (3) maintain
an adversarial journal-readiness verdict in `journal_readiness.md`;
(4) pursue promising directions with real experiments (E017+). Spencer
approved free editing of the paper, including folding in E014/E015/E016.
Target paper length: roughly 10 pages (was 14 before fold-in).*

## Iteration 1 (2026-07-22 18:26)

- Verified E014/E016 headline numbers against their READMEs (rho 0.76/0.63,
  -0.02; angles 71-88 deg; oracle 0.99) before pasting. E015 numbers already
  verified in `v2_density_law.md`.
- Folded all three into `qce2027_paper.tex`: E015 lemmas + conditioning
  remark after Corollary 1; E014 paragraph + `fig:argmax` in section 5.4;
  E016 replaces the Discussion's chicken-and-egg sentence; experiment range
  now E001--E016. Copied `argmax_vs_fidelity.png` into `Figures/generated/`.
- Compiles clean; now 15 pages. Committed in paper repo (`e09a5dc`).
- Updated STATUS.md: fold-in decision resolved, convention changed.

Next: full error-hunt pass (verify every number in the paper against
results.csv files) and begin trimming toward 10 pages; then explainer pass.

## Iteration 2 (2026-07-22 18:32)

- Ran a three-agent number audit of the whole paper against the experiment
  records (every claim, every table cell, plus re-derivation of all proofs).
- All mathematics re-derived correctly, including the new V2 lemmas. Most
  numbers confirmed exactly (regret table byte-identical to E004's, slack
  5.8x max, ratio 1.505, all timing entries, all fitted-shape stats).
- One real error found and fixed: "slice sums up to ~1e8" was wrong twice
  over (E012 truth: mean of per-instance maxima 9.2e8, absolute max 2.3e10;
  STATUS's "9e8" was a mean mislabeled as an "up to").
- Applied 15 precision/scoping fixes (see paper commit `b47de04`), including
  adding the two missing lemma proofs and fixing a notation collision
  (s_2 vs spin s_j, spins now sigma_j).
- Paper compiles clean; 16 pages now (proofs added). Trimming to ~10 pages
  is the next major chunk, likely: cut the full-page glossary table, merge
  figures, tighten Sec. 5 prose.

Next: length trim toward 10 pages, then explainer rewrite pass.

## Iteration 3 (2026-07-22 18:50)

- Trimmed the paper 16 -> 12 pages without cutting results: glossary table
  halved (standard QAOA notation folded into its caption), whyhom and
  trajectory-PCA paragraphs tightened, intro question paragraph condensed,
  timing paragraph condensed, 10pt/1in layout, figures at 0.8 linewidth,
  small bibliography font. Paper repo commit follows iteration 2's audit.
- Judgment call recorded: going below ~12 single-column pages means cutting
  actual results (candidates: drop the glossary entirely, drop fig:depth or
  fig:bound, compress section 5.4). In the target IEEEtran two-column
  format the current content is roughly 9-10 pages, so the 10-page goal is
  effectively met under the venue format; parked as a Spencer decision.

Next: explainer completeness pass (goal 2), then journal_readiness.md.

## Iteration 4 (2026-07-22 19:00)

- Explainer updated to full rewrite-from-scratch completeness: added the proxy
  recursion and expectation formulas (section 1), new section 4e (the V2
  codegree lemmas and the triangle-conditioning mechanism, plain language),
  the E014 argmax-vs-fidelity Q&A, rewrote section 7 (E016 closes the cheap
  frame door; V2 open piece narrowed; polynomial-N estimator listed), new
  section 9 (the paper's section-by-section skeleton with the related-work
  positioning map) and section 10 (key-numbers table with experiment
  provenance, post-audit values), cheat-sheet rows for delta_i, T, A_jk, tau,
  and post-audit number corrections throughout (slack median 4x/max 5.8x,
  regret 0.028-0.055).

Next: journal_readiness.md (goal 3), via adversarial referee agents.

## Iteration 5 (2026-07-22 19:10)

- Ran three adversarial referee reviews (novelty/significance, statistics/
  methodology, claims/presentation). Consensus: MAJOR REVISION at IEEE
  TQE/QST tier; no referee expected conclusions to reverse. Verdict and
  ranked gap list written to research/journal_readiness.md.
- Top gaps: (1) family-demeaned re-analysis + bootstrap CIs for all pooled
  Spearman rhos (E017, planned next); (2) grid-ceiling validation vs.
  continuous optimization + a literature fixed-angle baseline (E018);
  (3) any theory bound linking leakage to regret (hard, highest impact);
  (4) demonstrated (not extrapolated) crossover at n=24-28 (needs HPC GPU,
  parked); (5) exact Var(E[T|c]) for the density law.
- Applied all mechanical claim-discipline fixes the presentation referee
  found (paper commit 3e767d7), including scoping "never vacuous" (Fig. 2's
  extreme-ramp points exceed the trivial bound), the abstract's norm claim,
  and the 9-vs-6 model count.

Next: E017 (bootstrap CIs + within-family correlations, analysis-only,
local CPU) to close the biggest statistics gap.

## Iteration 6 (2026-07-22 19:30)

- E017 (statistics hardening) run and committed: bootstrap CIs, cluster
  bootstrap, within-cell and demeaned Spearman for E012/E014; CI-carrying
  7-point ranking rhos (E004); per-instance depth ratio 1.54 +- 0.26 (E011
  vs E003, 70 matched runs); 26/28 sign test p = 1.5e-6 (E010).
- Substantive finding: the fidelity-regret correlations were family-
  confounded in both directions (p=1 pooled 0.39 is mostly between-family,
  within-cell 0.07; p=3 pooled -0.02 hides within-cell 0.46). Paper's E014
  paragraph and figure caption rewritten honestly; argmax-displacement claim
  survives every control and is now the paper's most robust statistic.
- Paper commit 8418bd9; experiment range now E001-E017; 13 pp.
- journal_readiness gaps 1 and 6 marked done.

Next: E018 (grid-ceiling validation vs. continuous optimization + a
published fixed-angle baseline), local CPU at n=12-14.
