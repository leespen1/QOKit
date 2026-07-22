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

## Iteration 7 (2026-07-22 19:45)

- E018 (ceiling validation + transfer baseline) run on 140 instances:
  ceilings are tight (p=1 gap <= 0.0007; p=3 ramp gap ~0.002; the linear-
  ramp restriction itself costs ~0.007 mean, up to 0.05 sparse). The new
  parameter-transfer baseline (mean true argmax of ten brute-forced
  ER(0.5) n=12 instances, applied verbatim everywhere) BEATS the exact
  compression on 134/140 instances at p=1 and 140/140 at p=3 (pooled
  transfer regret ~0.014 / ~0.008 vs proxy 0.032 / 0.083). The paper now
  validates its ceilings in the setup, reports the baseline in section 5.3,
  and delimits the proxy's niche in the abstract. Honest but significant
  reframing: on concentrated random ensembles, per-instance parameter
  setting adds little over transfer.
- Literature agent (gap 7): no peer-reviewed instance of entrywise
  N-fitting exists; the documented origin is the collaboration's own public
  G-RIPS 2024 report, now cited (with Khairy 2020 / Shaffer 2023 for the
  general surrogate-fit pattern). Table 3 framing rewritten.
- Paper commit fb036e1; experiment range E001-E018; 13 pp; compiles clean.
- journal_readiness gaps 2 and 7 marked done. Remaining open: theory-regret
  bound (gap 3), HPC crossover demo (gap 4, parked), Var(E[T|c]) (gap 5),
  vector figures (gap 8).

Next: attempt gap 5 (exact conditioning correction Var(E[T|c]) for G(n,p)),
the most tractable remaining science; if it stalls, do gap 8 (vector
figures).

## Iteration 8 (2026-07-22 ~20:00)

- E019: closed the density-law open problem in substance (readiness gap 5).
  Derived and machine-verified four exact moment identities; the key one is
  Cov(T, S^2) = Var(T) (a quadratic channel the crude 36 tau^2/m estimate
  missed, alive even in triangle-free graphs). The resulting L2-projection
  bound V2 <= Var(T) - v'G^{-1}v is rigorous, polynomial-time, and captures
  91-100% (mean 97%) of the exact conditioning correction over 140
  instances; ER asymptotic V2/Var(T) -> 1 - p(1+4p-2p^2)/(1+6p^2-4p^3)
  (0.375 at p=1/2 vs 0.352 measured). All identities asserted to 1e-9.
- Paper: new Proposition (quadratic conditioning bound) with proof in
  section 4; Discussion open problem narrowed; abstract upgraded to
  "derive to within a few percent". Commit 24a1c7f; 13 pp; compiles clean.

Next: digest E019 into STATUS/explainer, then remaining gaps: theory-regret
bound (gap 3, hard) or vector figures (gap 8, mechanical).

## Iteration 9 (2026-07-22 ~20:15)

- Digested E019 into STATUS.md and the explainer (new density-law entry,
  key-numbers rows for E017/E018/E019).
- Figure-production pass (readiness gap 8) done via subagent: all six paper
  figures re-emitted as vector PDFs at print font sizes; "gate criterion"
  annotation removed from the ranking figure; "Normal" renamed "Gaussian" in
  the fitted-shape figure. QOKit commit cb12e3cf, paper commit 44b8dd7;
  paper still 13 pp, compiles clean.
- Readiness gaps now: 3 (theory bound linking leakage to regret; hard) and
  4 (HPC crossover demo; parked for Spencer) remain open; 1, 2, 5, 6, 7, 8
  closed.

Next: attempt gap 3 (a p=1 landscape-error -> regret bound through peak
flatness), the last open science this loop can reach; then a refreshed
referee verdict near the end of the 24 h window.

## Iteration 10 (2026-07-22 ~20:45)

- E020: closed the last open science gap this loop can reach (gap 3). Proved
  a two-point regret certificate connecting the leakage calculus to regret
  (regret <= eps(theta*) + eps(theta-hat), eps controlled by leakage and
  cost variances), machine-verified on 280 instance-depth pairs, honest
  tightness median 10-11x; the looseness is the quantitative form of the
  argmax-transfer thesis. New Proposition + proof in paper section 4.
- Bonus practical finding: the certificate's normalized proxy objective
  changes the selected grid point on 140/140 instances and HALVES depth
  regret (p=3: 0.080 -> 0.044, better on 134/140) at zero cost; p=1 mildly
  prefers the unnormalized convention. New paragraph + bolded rule in
  section 5.4; abstract updated. Paper commit 514b221; 14 pp.
- All eight journal_readiness gaps now closed or parked (gap 4 needs HPC).

Next: digest E020 into STATUS/explainer; then re-run the referee panel for
a refreshed verdict, and a final consistency pass (page count crept to 14;
one more trim).

## Iteration 11 (2026-07-22 ~21:15)

- Refreshed referee review returned: verdict upgraded to MINOR REVISION.
  All three original objections judged addressed (the reframe from method
  to anatomy accepted; certificate closes the wrong-error objection; the
  statistics now "referee-proof in structure").
- Applied its five desk edits (paper commit follows): abstract rewritten to
  ~200 words; "How to set parameters" decision rule placed beside the
  recipe; normalization convention pinned in setup, Table 2 caption, and a
  certificate remark; "So when does it work?" boxed answer opens the
  Discussion; glossary cell + wording hygiene.
- journal_readiness.md verdict updated; only Spencer-owned items block
  submission (author list, acknowledgments, frozen archive/DOI, venue).

Next: sweep the explainer for consistency with today's paper changes
(abstract, decision rule), then final wrap-up passes within the window.

## Iteration 12 (2026-07-22 ~21:45)

- E021: proved-by-verification that the whole section-4 theory extends
  verbatim to integer-weighted MaxCut (weighted codegrees, triangle and
  4-cycle weight products): all identities to 1e-9 on 140 weighted
  instances, quadratic bound capture 93-100% (mean 0.975), cubic law to
  0.4%. Also established the structural boundary: continuous weights
  collapse every cost class to a bitstring-complement pair, so the
  compression itself becomes vacuous (checked on all 7 families).
- Paper: new weighted-MaxCut remark in section 4; Limits rescoped
  (experiments unweighted, theory weighted-ready). Commit 953d914; 14 pp.
- STATUS and explainer digested.

Next: explainer consistency sweep vs. today's paper changes (abstract
rewrite, decision rule, certificate), then final wrap-up.

## Iteration 13 (2026-07-22 ~22:45)

- Fresh two-agent verification of everything added today: all new math
  re-derives correctly (including the exact identity 1-||phi||^2 = sum of
  lambda^2 for the unnormalized certificate term). Eight paper fixes
  applied (commit 1b6358d): grips2024 collaborator initials corrected
  (S. Kakuta, K. Sakurai); weighted-remark replacement dictionary made
  precise (second moments map to sum w^2 / sum w^4, not W); contributions
  list synced with the certificate, conditioning bound, transfer baseline,
  and normalization rule; transfer p=3 regret restated as 0.006 under the
  stated ceiling convention; capture/residual ranges harmonized.
- Explainer fully synced with today's paper (QOKit commit 723dc139):
  new sections for the conditioning bound, weighted extension, and regret
  certificate with proof ideas; corrected E014/E017 numbers; skeleton and
  decision-rule/"So when does it work?" content; abstract-shape note.

Next: overnight cadence. Remaining productive options: attempt sharpening
the certificate (correlated-error analysis, E022); a final full-paper
re-audit near the end of the window; keep STATUS decisions current.

## Iteration 14 (2026-07-22 ~23:30)

- E022: exact signed decomposition regret = Delta-e - margin (identity
  asserted 280/280). The certificate's 10x slack has three named parts:
  uniform overprediction by the normalized proxy (e<0 at 277/280 argmax
  points, so signed errors partially cancel, median factor 0.6), a winner's
  curse concentrating error at the proxy's own argmax (4-5x the error at
  the true optimum), and the proxy's internal margin absorbing a third of
  the residual tilt. Paper's certificate paragraph now states this
  (commit 430027f); experiment range E001-E022.
- Sharp insight for future theory: regret bounds must model the selection
  effect, not pointwise error. This is the mechanistic reason no pointwise
  norm predicts regret (E012/E014).

Next: overnight lower-intensity cadence; near end of window, final
full-paper audit + wrap-up summary for Spencer.

## Iteration 15 (2026-07-23 ~00:15)

- Housekeeping: proposed_paper_additions.md marked superseded;
  theory_compression.tex flagged out-of-sync (Spencer decision 5);
  STATUS decision 3 (headline framing) marked resolved-in-effect.
- E023 (negative): the E022 overprediction is not predictable from norm
  loss (per-instance rho -0.26 to 0.60) and a pooled linear correction
  worsens regret on every moved instance (0/126 better, both depths).
  Third confirmation of "never calibrate values to fix argmaxes"; the
  winner's-curse remainder is a selection effect. One sentence added to
  the paper's normalization paragraph (commit dab97cd); range E001-E023.
- Page count 15 after all additions; the venue-format argument (two-column
  IEEE ~10-11 pp) still applies; further single-column trimming would cut
  results, parked as Spencer's call.

Next: overnight; near end of window, final full re-audit + wrap-up
summary. Science backlog now: correlated-error/selection-aware regret
theory (hard), HPC crossover demo (parked).
