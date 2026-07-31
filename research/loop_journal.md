# 24-hour improvement loop journal

*Loop started 2026-07-22 18:26 JST; hard stop after 2026-07-23 18:26 JST.
Spencer's directives: (1) hunt and fix paper errors; (2) make
`paper_explainer.md` complete enough to rewrite the paper from; (3) maintain
an adversarial journal-readiness verdict in `journal_readiness.md`;
(4) pursue promising directions with real experiments (E033+). Spencer
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
  Spearman rhos (E033, planned next); (2) grid-ceiling validation vs.
  continuous optimization + a literature fixed-angle baseline (E034);
  (3) any theory bound linking leakage to regret (hard, highest impact);
  (4) demonstrated (not extrapolated) crossover at n=24-28 (needs HPC GPU,
  parked); (5) exact Var(E[T|c]) for the density law.
- Applied all mechanical claim-discipline fixes the presentation referee
  found (paper commit 3e767d7), including scoping "never vacuous" (Fig. 2's
  extreme-ramp points exceed the trivial bound), the abstract's norm claim,
  and the 9-vs-6 model count.

Next: E033 (bootstrap CIs + within-family correlations, analysis-only,
local CPU) to close the biggest statistics gap.

## Iteration 6 (2026-07-22 19:30)

- E033 (statistics hardening) run and committed: bootstrap CIs, cluster
  bootstrap, within-cell and demeaned Spearman for E012/E014; CI-carrying
  7-point ranking rhos (E004); per-instance depth ratio 1.54 +- 0.26 (E011
  vs E003, 70 matched runs); 26/28 sign test p = 1.5e-6 (E010).
- Substantive finding: the fidelity-regret correlations were family-
  confounded in both directions (p=1 pooled 0.39 is mostly between-family,
  within-cell 0.07; p=3 pooled -0.02 hides within-cell 0.46). Paper's E014
  paragraph and figure caption rewritten honestly; argmax-displacement claim
  survives every control and is now the paper's most robust statistic.
- Paper commit 8418bd9; experiment range now E001-E033; 13 pp.
- journal_readiness gaps 1 and 6 marked done.

Next: E034 (grid-ceiling validation vs. continuous optimization + a
published fixed-angle baseline), local CPU at n=12-14.

## Iteration 7 (2026-07-22 19:45)

- E034 (ceiling validation + transfer baseline) run on 140 instances:
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
- Paper commit fb036e1; experiment range E001-E034; 13 pp; compiles clean.
- journal_readiness gaps 2 and 7 marked done. Remaining open: theory-regret
  bound (gap 3), HPC crossover demo (gap 4, parked), Var(E[T|c]) (gap 5),
  vector figures (gap 8).

Next: attempt gap 5 (exact conditioning correction Var(E[T|c]) for G(n,p)),
the most tractable remaining science; if it stalls, do gap 8 (vector
figures).

## Iteration 8 (2026-07-22 ~20:00)

- E035: closed the density-law open problem in substance (readiness gap 5).
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

Next: digest E035 into STATUS/explainer, then remaining gaps: theory-regret
bound (gap 3, hard) or vector figures (gap 8, mechanical).

## Iteration 9 (2026-07-22 ~20:15)

- Digested E035 into STATUS.md and the explainer (new density-law entry,
  key-numbers rows for E033/E034/E035).
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

- E036: closed the last open science gap this loop can reach (gap 3). Proved
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

Next: digest E036 into STATUS/explainer; then re-run the referee panel for
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

- E037: proved-by-verification that the whole section-4 theory extends
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
  certificate with proof ideas; corrected E014/E033 numbers; skeleton and
  decision-rule/"So when does it work?" content; abstract-shape note.

Next: overnight cadence. Remaining productive options: attempt sharpening
the certificate (correlated-error analysis, E038); a final full-paper
re-audit near the end of the window; keep STATUS decisions current.

## Iteration 14 (2026-07-22 ~23:30)

- E038: exact signed decomposition regret = Delta-e - margin (identity
  asserted 280/280). The certificate's 10x slack has three named parts:
  uniform overprediction by the normalized proxy (e<0 at 277/280 argmax
  points, so signed errors partially cancel, median factor 0.6), a winner's
  curse concentrating error at the proxy's own argmax (4-5x the error at
  the true optimum), and the proxy's internal margin absorbing a third of
  the residual tilt. Paper's certificate paragraph now states this
  (commit 430027f); experiment range E001-E038.
- Sharp insight for future theory: regret bounds must model the selection
  effect, not pointwise error. This is the mechanistic reason no pointwise
  norm predicts regret (E012/E014).

Next: overnight lower-intensity cadence; near end of window, final
full-paper audit + wrap-up summary for Spencer.

## Iteration 15 (2026-07-23 ~00:15)

- Housekeeping: proposed_paper_additions.md marked superseded;
  theory_compression.tex flagged out-of-sync (Spencer decision 5);
  STATUS decision 3 (headline framing) marked resolved-in-effect.
- E039 (negative): the E038 overprediction is not predictable from norm
  loss (per-instance rho -0.26 to 0.60) and a pooled linear correction
  worsens regret on every moved instance (0/126 better, both depths).
  Third confirmation of "never calibrate values to fix argmaxes"; the
  winner's-curse remainder is a selection effect. One sentence added to
  the paper's normalization paragraph (commit dab97cd); range E001-E039.
- Page count 15 after all additions; the venue-format argument (two-column
  IEEE ~10-11 pp) still applies; further single-column trimming would cut
  results, parked as Spencer's call.

Next: overnight; near end of window, final full re-audit + wrap-up
summary. Science backlog now: correlated-error/selection-aware regret
theory (hard), HPC crossover demo (parked).

## Iteration 16 (2026-07-23 ~01:30)

- E040 (n=16 stability, 35 instances): all three parameter-space findings
  replicate one size up (displacement-regret rho 0.84/0.76; normalization
  rule 0.090->0.060 at p=3, better 32/35, p=1 keeps unnormalized edge;
  signed anatomy 69/70 overprediction, winner's curse 4-9x). Two stability
  sentences added to the paper (commit 24a1e1b); range E001-E040.

Next: overnight idle stretch; final comprehensive re-audit and Spencer
wrap-up digest in the last hours of the window.

## Iteration 17 (2026-07-22 ~22:25)

- Repo health check: full JuliaQAOA test suite passes (exit 0, all test
  summaries Pass) after the day's additions. No source changes were made by
  the loop (experiments are self-contained), as expected.
- Entering overnight idle cadence: hourly heartbeats; the final
  comprehensive audit and the Spencer wrap-up digest are scheduled for the
  last hours of the window (~14:00-18:00 JST on 2026-07-23).

## Iteration 18 (2026-07-22 ~23:00)

- Explainer coherence pass: the ladder subsection relabeled (it follows the
  certificate section, matching the paper's order); the certificate's "how
  to read it" now carries the E038 signed anatomy, the E039 negative, and
  the E040 n=16 stability in prose (they were table-only).
- Overnight: idle heartbeats; final audit + Spencer digest in the last
  hours of the window.

# Loop 2 (72-hour weekend run)

*Loop started 2026-07-31 19:12 JST; hard stop 2026-08-03 ~19:12 JST (Spencer
extended the originally-requested 24 h to a full weekend mid-launch). Same
four directives as loop 1: (1) hunt and fix paper errors; (2) explainer
complete enough to rewrite the paper from; (3) honest journal-readiness
verdict, now judged explicitly against a low-to-mid journal tier;
(4) pursue promising directions with real experiments (E041+). No trivial
questions to Spencer. Plan: ~/.claude/plans/i-want-you-to-zippy-seal.md.*

## Iteration 1 (2026-07-31 19:12)

- Housekeeping: adopted the orphaned E041 (recipe composition: sampled N x
  normalized objective) left by loop 1 — its full run died at 123/140
  instances at the window boundary. README written (Answer: PENDING), full
  run relaunched in background (~15 min). Stray E012/E014/E016 smoke CSVs
  committed (QOKit commit 0bedf77d); tree clean.
- Length question SETTLED: a scratch IEEEtran conference-class build of the
  unmodified paper compiles clean at 9 pp two-column. The 10 pp target is met
  with no cuts; 15 pp single-column is a formatting artifact. STATUS updated
  (decision 2 and the headline paragraph).
- In flight at iteration close: E041 full run (background julia, ~15 min);
  two audit agents (numbers-vs-records and math re-derivation), prioritizing
  the post-July-22 additions (E033-E040 content, certificate, conditioning
  bound, weighted remark).

## Iteration 2 (2026-07-31 19:35)

- E041 complete (280 rows): the recipe composes and MORE — sampled N (S=10)
  + normalized is the best of all four {exact,sampled}x{raw,norm} cells at
  BOTH depths, beating exact+normalized (p=1 mean regret 0.021 vs 0.047,
  sign p=8e-29; p=3 0.024 vs 0.044, p=2e-15), uniformly across all 7
  families. Exact-N cells reproduce E036/E040 (p=3 0.080->0.044; p=1 exact
  prefers raw). README digested; claim held OUT of the paper pending
  mechanism.
- E042 launched (mechanism: noise vs estimator bias): R=10 replicates +
  replicate-averaged N + S-sweep {3,10,30,100,300}, n=12, same seeds.
  Smoke passed; full run in background.
- Audit agents (numbers, math) still running.

## Iteration 3 (2026-07-31 ~20:20)

- Dual audit returned and applied: ~79 claim groups checked against records,
  68 confirmed, 11 mismatches fixed (one substantive: norm vs squared-norm
  7.5x; the rest precision/labeling: SE range 0.003-0.008, winner's curse
  3.8-4.7x, margin ~40%, 557/560 denominators, fitparadox caption
  composition, n=16 correlation convention disclosed, transfer claim scoped
  to proxy-only methods). Math audit: NO broken theorems, all identities
  machine-verified; fixed an undefined epsilon, a missing 1/c_opt unit in
  the certificate addendum, added the O(beta^2 gamma) cancellation sentence
  to Cor. 6, singular-G qualifier and bound-not-ratio asymptotic in Prop 10.
  Local paper commit ea8f56c (compiles clean, 15 pp single-column).
  Deferred: T/S/v/p/G symbol collisions (dedicated notation pass later).
- E042+E043 digested and committed: the sampled-N gain is NOISE acting
  through the normalized objective (E042: averaged-N regresses to exact,
  S=3 best, raw indifferent); at p=1 generic iid noise suffices, at p=3
  only the sampler's structured (marginal-preserving) noise helps (E043).
  Strong practical + conceptual finding, held out of the paper until the
  reconciliation below lands.
- **DISCOVERY — the paper repo has two diverged lines.** Remote
  ClaudeResearch is 37 ahead: the July editorial line (E017-E032 in ITS
  numbering = 017_error-directions..032_transfer-scale-audit; QCE 9 pp cut
  d161076; QINP 17 pp journal manuscript de6bbaf/2953996; Spencer-attributed
  decisions: Sud reconciliation, weighted-arc reframe, E032 audit). Local is
  30 ahead: loop-1's line (June-13 base, never pulled; E014-E024 fold-in in
  loop-1 numbering = 017_statistics-hardening..024_n16-stability;
  certificate, conditioning bound; today's audit fixes). Experiment numbers
  017-027 are DOUBLED in research/experiments/ with distinct slugs; E-refs
  ambiguous from 017 up. Push to the paper repo is blocked (no force-push).

**Reconciliation plan (next iteration's main task):**
1. In QOKit: renumber the loop-1-derived series 017_statistics-hardening ->
   033, 018_ceiling-validation -> 034, 019_conditioning-correction -> 035,
   020_regret-certificate -> 036, 021_weighted-maxcut -> 037,
   022_signed-regret-decomposition -> 038, 023_bias-corrected-objective ->
   039, 024_n16-stability -> 040, 025_sampled-plus-normalized -> 041,
   026_sampled-N-mechanism -> 042, 027_synthetic-noise-control -> 043
   (git mv; July series keeps its numbers since the REMOTE paper cites
   them). Update refs in loop_journal, paper_explainer, STATUS,
   journal_readiness, and the loop-1 experiment READMEs/scripts.
2. In the paper repo: treat remote (QINP head 2953996) as the canonical
   editorial line. Port loop-1's + today's science (E014-E016 fold-in,
   statistics hardening, certificate + signed anatomy, conditioning bound,
   weighted remark, n=16 stability, audit fixes) INTO it as new commits
   with the renumbered E-refs, using the local line as source material.
   Merge strategy: merge -s ours of the local branch afterward to join
   histories without discarding either (local commits stay reachable).
3. Correct STATUS: my 2026-07-31 "15 pp / 9 pp two-column" headline
   described the stale local line; remote already has a 9 pp QCE cut AND a
   17 pp QINP manuscript. The 10 pp target discussion must name the actual
   files.
