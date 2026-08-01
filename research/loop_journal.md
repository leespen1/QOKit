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

## Iteration 4 (2026-07-31 ~20:50)

- Renumber executed: loop-1 series git-mv'd to 033-043; references updated
  in loop_journal (pre-iteration-3 bulk + hand-fixed mixed passages),
  paper_explainer, journal_readiness, STATUS (loop-1 digest blocks only;
  July blocks untouched), and the renumbered dirs' own files. A glob
  mistake briefly touched July dirs 031/032 — caught and reverted before
  commit. STATUS gains an "E-number collision resolved" section with the
  full map. QOKit pushed.
- Port agent launched in the paper repo: branch reconcile-loop1 from
  origin/ClaudeResearch (QINP manuscript = editorial flagship), porting
  loop-1's certificate/conditioning/normalization/statistics/E014-E016
  science with renumbered refs, honoring the remote line's owner-approved
  framings (weighted reframe, Sud tone). No pushes; I review before merge.

## Iteration 5 (2026-07-31 ~21:05)

- Paper reconciliation COMPLETE and pushed (2953996..1e72f46, fast-forward;
  no force-push; loop-1 commits joined via ours-merge, reachable forever).
  Four port commits put loop-1's science into qinp_paper.tex (17->20 pp,
  clean builds): dual-audit fixes to shared claims; conditioning bound as
  the primary V2 statement (remote's quartic-split compressed to a pointer;
  capture set to the README-backed 91-100%/97%); certificate + signed
  anatomy + normalization rule (with a bridging sentence reconciling
  "normalize at depth" with the remote's "normalization is an instrument"
  finding); statistics hardening + ceiling validation + transfer
  replication merged into the remote's own framing. Weighted remark placed
  in theory, rescoped to MOTIVATE the binned-proxy arc; a genuine remote
  imprecision fixed (unweighted identity DOES survive integer weights).
- Key scoping kept from the remote line (harsher than loop-1's story):
  at p=10-20 every proxy variant collapses (0.07-0.37) while ramp transfer
  holds <=0.03; "sampled N" is best-proxy-at-p<=3, not a universal recipe.
  E041-E043's noise findings must be read within that scope when they
  enter the paper.
- STATUS "Where we are" rewritten: two manuscripts named (QCE 9 pp — 10 pp
  target met; QINP 20 pp journal flagship), reconciliation recorded.

## Iteration 6 (2026-07-31 ~21:35)

- Referee panel returned: **MINOR REVISION equivalent at QINP tier; the
  science is submission-ready, the package is not.** Direct answer to
  Spencer's suitability question: YES for a low-to-mid journal, QINP a good
  fit, no further research required; E041-E043 is the strongest optional
  upgrade (with a consistency obligation: the paper's "normalize at depth"
  recipe is dominated by sampled+normalized per E041, whose quarantine
  condition is now met by E042). journal_readiness.md rewritten (commit
  7074f588). Top gaps: (1) declarations+Zenodo (Spencer-owned), (2) the
  E041-E043 fold-in decision, (3) vector figures, (4) addendum machine
  check + n=16 sentence, (5) front matter.
- Gap 4 closed: verify_addendum.jl machine-checks the unnormalized
  addendum at 714 points (norm identity, decomposition, full bound, 1e-10;
  worst utilization 0.85). E036 README updated; paper provenance tag now
  cites the check. Both repos pushed.
- Gap 3 delegated: vector-figure agent running on the QINP build.
- Explainer rebuild agent still running.

Next: when the explainer lands, decide the E041-E043 fold-in (gap 2) —
one scoped paragraph in section 6.1 + a Discussion pointer, keeping the
remote line's "best proxy variant at p<=3" scoping.

## Iteration 7 (2026-07-31 ~21:55)

- Explainer rebuild landed (commit 0e8fa6b6 + flag fix; pushed): tracks the
  QINP manuscript end to end (July arcs incl. the harsh-depth E019 result,
  binned weighted proxy, E032 retraction; ported loop-1 science; section 8
  = E041-E043 frontier marked "not yet in the paper"; ~50-row key-numbers
  table with provenance; 10,348 words).
- Stale-records cleanup: seven READMEs aligned with CSV-backed audit
  numbers (034 transfer 0.006; 014 rhos 0.74/0.65; 011 +11-33%/25-57%;
  012 slice-sum mean-vs-max; 035 capture min 0.9046; 038 curse 3.8-4.7x;
  004 squared norm).
- Noted for the next paper edit: QINP says capture "91-100%"; CSV min is
  0.9046, so it should read 90-100% — fold into the E041-E043 edit.
- Figure agent still running in the paper repo (holding off on paper edits
  until it lands to avoid collisions).

## Iteration 8 (2026-07-31 ~22:15)

- Readiness gap 2 closed: E041-E043 folded into qinp_paper.tex as one
  scoped paragraph ending Sec. 6.1 (composition beats every combination at
  both depths; mechanism = in-set noise disrupting the winner's curse;
  scoped best-proxy-at-p<=3, transfer still wins in-family). Capture floor
  corrected to 90% while in there. Gap 3 closed by the figure agent (all
  six QINP figures vector, commit c2857c3). Gap 4 closed last iteration.
  journal_readiness updated: only gap 1 (Spencer: declarations + Zenodo)
  and gap 5 (front matter) remain.
- Explainer section 8 and key-numbers rows updated to reflect the fold-in.

Next: gap 5 front-matter pass (abstract compression, keywords, watermark,
bibliography re-verification), then fresh science (E044 candidates: S and
n scaling of the noise-regularization effect; selection-aware certificate
sharpening using the E038 anatomy).

## Iteration 9 (2026-07-31 ~22:40)

- Front-matter agent launched on qinp_paper.tex (abstract compression,
  keywords, watermark sweep, bibliography re-verification, Springer
  declarations block with TBD-Spencer placeholders).
- E044 launched (n=16 replication of the composed-recipe headline, plus
  the S=3 cell): smoke passed, full run in background (~40 min).

## Iteration 10 (2026-07-31 ~23:00)

- Front-matter pass landed and pushed (paper commit e429241): abstract
  215->197 words (honest framing intact), Springer keywords + MSC block,
  draft date removed, declarations section with TBD-Spencer placeholders,
  montanezbarrera2025 title fixed; 21 pp, clean build, zero undefined refs.
  All 35 citations resolve; every bibitem cited.
- Ten flagged references (incl. a suspect lotshaw2023 venue) now being
  verified against arXiv/publisher records by a web-verification agent.
- E044 (n=16 composed-recipe replication) still running.

## Iteration 11 (2026-07-31 ~23:25)

- Bibliography fully verified: all 11 flagged entries checked against
  publisher/arXiv records (paper commit 76282fe, pushed). One fix
  (montanezbarrera2025 published title); lotshaw2023 suspicion unfounded
  (PRA 108, 042411 is correct). Every reference in the manuscript is now
  record-verified. All non-Spencer submission items are done.
- E044 still in prep (serial exact-N at n=16 is the bottleneck); expected
  to land within the hour.

## Iteration 12 (2026-07-31 ~23:50)

- E044 landed: the composed recipe replicates at n=16, 34/35 at both
  depths (p=1e-9), margin WIDENS with n (0.018/0.024 vs exact-normalized
  0.061/0.060); S=3 strongest (0.0067 at p=1). Exact-norm cell reproduces
  E040. Digested into README, one clause in qinp Sec. 6.1 (% E044,
  pushed), and the explainer key-numbers table.
- The E041-E044 arc is now a complete, size-stable, mechanism-explained
  result fully inside the paper.

Next: E045 (selection-aware certificate sharpening via the E038 signed
anatomy: can regret <= max(0, Delta-e) - margin with Delta-e estimated
from the norm-visible part give a tighter certified bound?), and an
overnight cadence after that.

## Iteration 13 (2026-08-01 ~00:30)

- E045 agent derived the margin-subtracted certificate and committed
  run.jl + smoke; the full run is still executing (trajectory computation
  is the heavy part). Polling until results.csv lands, then README +
  digest. Overnight cadence.

## Iteration 14 (2026-08-01 ~01:40)

- E045 landed and digested: the margin-subtracted certificate is rigorous
  and never worse (140/140) but closes little (median 10.2->9.5 / 11.3->
  10.6); margin ~0.03 vs eps sum ~0.46; error cosines (-0.43/-0.05) show
  limited headroom for two-point refinements. Practitioner-uniform
  no-oracle variant validated as corollary. Honest negative, equal
  billing: README + one paper sentence (% E045, pushed). The selection-
  effect modeling gap stays open and is now sharply delimited.
- Overnight: next big item is tomorrow's full-paper coherence audit; the
  explainer needs an E045 row at the same time.

## Iteration 15 (2026-08-01, morning)

- Spencer directive: spend remaining subscription tokens on the work.
  Approved: estimator arc + verification gauntlet, workflows OK.
- Launching: E046 design/prototype agent (doubly-sampled poly-time N
  estimator: MCMC class members x sampled neighbor profiles; the E042/E043
  noise findings suggest crude estimation may WIN, not just suffice);
  and a reproduction-sweep workflow over all experiment dirs (rerun
  locally feasible ones from committed scripts, diff results; analysis-
  only verification for HPC-scale ones). Referee gauntlet queued after
  the sweep so its fixes land first.

## Iteration 16 (2026-08-01 ~05:40)

- LESSON: the 45-agent reproduction workflow + E046 agent hit the session
  limit within a minute; all 46 agents died with zero output (~235k tokens
  wasted). Spencer flagged it. Strategy revised: NO wide fan-outs;
  sequential small agents; inline main-loop work preferred (far cheaper).
  Limit resets 12:10am JST per the error.
- Inline sweep started: E001 REPRODUCED_FULL (identities regenerate,
  settling audit U3), E015 REPRODUCED_FULL (byte-identical results.csv).
  Tracking in research/reproduction_sweep_2026-08-01.md.
- E046 estimator agent died pre-output; will relaunch as a single agent
  when limits allow (probe first with a tiny agent).

## Iteration 17 (2026-08-01 ~06:50)

- Limits back (probe ok). E046 estimator relaunched as a single agent.
- Sweep: E016 byte-identical (210 jobs). E033-E040 marked ANALYSIS_VERIFIED
  on the strength of the 2026-07-31 audit's CSV recomputations (036 also
  had its addendum check rerun this loop). Ledger at 19/45.

## Iteration 18 (2026-08-01 ~08:15)

- E046 preview (run still in its n=20 scaling phase): the poly-time
  estimator VALIDATES — p=1 norm regret 0.011 (best of all variants;
  samp10 0.019, exact 0.047); p=3 0.030 vs samp10 0.025, exact 0.045.
  Within the pre-set ~0.01 AR criterion: the paper's open question is
  closing positively. Awaiting the agent's README/timing before any
  paper edit.
- Sweep: E004 (gate) byte-identical earlier this morning; E005+E006
  reruns still executing.

## Iteration 19 (2026-08-01 ~10:15)

- **E046 COMPLETE: the paper's open question closes positively.** A fully
  polynomial pipeline (Wang-Landau density of states -> flat-histogram
  MCMC class members -> K=200 random d-subsets per profile row ->
  estimated counts + normalized objective) matches the exponential
  enumeration sampler within 0.01 AR at n<=14 (p=1: BEST variant, mean
  0.011, better 49/70 vs samp10; p=3: 0.030 vs 0.025, within band),
  misses zero attained classes, and crosses below enumeration wall-clock
  at n=20 (0.96s vs 1.70s, 3 instances). QOKit commits 43de9c98/4755e46d.
  Paper updated in two places (contributions item, practical subsection)
  with prototype-grade scoping; 21 pp, clean build, pushed.
- The recipe's domain converts from a hard n<~30 ceiling to an
  extrapolation question. Raw objective collapses under estimated counts:
  the normalized objective is an essential partner (consistent with the
  whole E041-E044 mechanism story).

## Iteration 20 (2026-08-01 ~11:15)

- Four-referee gauntlet returned: 2 BLOCKERs (both my own mid-sentence
  %-comment bugs swallowing words in the rendered PDF), 7 MAJOR, 15 MINOR;
  11 fixes applied and pushed by the agent (incl. purging three leftover
  "only the analytical N survives" sentences that contradicted E046,
  softening "fully polynomial" to per-sweep-polynomial with a WL
  heuristic-convergence caveat + citation, aligning E011 numbers to the
  record, E-range pointers to E001-E046). One genuine number discrepancy
  found and fixed (E011 +13-30% -> +11-33%). Verdict unchanged: MINOR
  REVISION at QINP tier.
- Deferred items resolved by coordinator: "margin widens" softened to
  "at least the same margin" with the formal-test caveat; 49/70 sign
  count disclosed; abstract left as-is (stakes the calculus; the
  contributions list carries the new results). Remaining open: Spencer's
  declarations + Zenodo archive only.

## Iteration 21 (2026-08-01 ~12:30)

- **Reproduction sweep COMPLETE: 45/45, zero reproduction failures.**
  17 full local reruns all byte-identical (including E003's 840 p=20 runs
  and the E041-E046 series); the rest analysis-verified (every README
  headline recomputes from committed CSVs — the July 018-032 agent needed
  no HPC data at all). Three README-side imprecisions fixed, all
  conservative or tolerance-level, none in the paper. Summary paragraph
  in research/reproduction_sweep_2026-08-01.md. The evidence base
  regenerates.
- Verification-gauntlet arm of the token plan is now done (audits ->
  referee gauntlet -> reproduction sweep). Science arm done through E046.
  Remaining before the window closes: STATUS refresh, final coherence
  read of the paper as one document, wrap-up digest.

## Iteration 22 (2026-08-01 ~13:10)

- STATUS What-we-know refreshed with four new entries: the noise-
  regularization arc (E041-E044), the estimator closure (E046), the E045
  negative, and the reproduction sweep. Pushed.
- E047 launched (single agent): estimator hardening — (S,K)
  hyperparameter robustness at n=12, n=16 validation (does poly's p=1
  advantage keep growing with n?), and timing at n in {16,20,24}
  including Wang-Landau step counts.

## Iteration 23 (2026-08-01 ~14:40)

- E047 landed and digested inline (agent stood down cleanly after the
  wasteful wake-cycle was cut off). Verdict split by depth: at p=1 the
  estimator is robust and IMPROVES with n (n=16: 0.0084, better 29/35;
  K=50 best at n=12 with 0.0061 — noise mechanism again); at p=3 the
  deficit to the enumeration sampler is not repaired by any (S,K) and
  WIDENS at n=16 (+0.022, better only 7/35). n=24 estimation ~1.3 s,
  ~2.6M WL steps, zero missed classes at all sizes. Paper scoping
  corrected in both E046 sites: "lifting the ceiling at low depth",
  explicit depth-deficit numbers, "prefer the enumeration sampler at
  depth". Pushed both repos.
- Science threads now closed for the window. Remaining: explainer E047
  touch, final one-document coherence read, wrap-up digest.

## Iteration 24 (2026-08-01 ~15:30)

- Coherence read landed (paper commit bc32dfc, pushed): the paper reads as
  one document; stale E046-range/tags fixed, the n<=14-vs-n=16 validation
  contradiction reconciled (FLAG for Spencer: the text now says validated
  to n=16, matching its own numbers; E047's README frames n=16 as "one
  size beyond validation" — pick a framing), Background headings starred,
  cites tied. Numbers agree everywhere spot-checked.
- E048 launched (background julia, ~2-3 h): the n=18 point for the
  estimator's depth-split trends, closing E047's two-points caveat.

## Iteration 25 (2026-08-01 ~17:00)

- E048 landed and digested: at n=18 the estimator's p=1 advantage is a
  statistical tie (+0.001, 11w/10l) — "grows with n" was premature and is
  now corrected in the paper to "beats at n<=16, ties at n=18"; the depth
  deficit is confirmed on three sizes (0.005/0.02/0.026) and quoted as
  such. Validation frontier updated to n=18. Bonus: the composed recipe
  extends to n=18 (beats exact-normalized 20/21 at p=1, 17/19 decided at
  p=3), added to the paper's stability sentence. Zero missed classes
  again. Both repos pushed; 21 pp clean.
- All science and verification threads for the window are now closed.
  Remaining: long-idle cadence, then the Spencer wrap-up digest in the
  final hours (2026-08-03 afternoon JST).

## Iteration 26 (2026-08-01 ~18:20)

- QCE conference cut synced with the flagship (paper commit c4b8cb0,
  pushed): applicable audit fixes applied (7 present, 4 absent-from-cut
  verified), one compact E041-E048 paragraph added at conference altitude,
  four stale "estimator is open / exponential essential" statements
  purged, Limits + repro footnote updated. **Exactly 10 pp** — Spencer's
  target met with the new results included.
- Memories saved: paper-repo structure (two manuscripts + numbering
  split) and agent-budget discipline (no wide fan-outs; the 45-agent
  lesson).
- ALL planned work for the window is complete: both manuscripts current
  and consistent, explainer current through E048, readiness at minor
  revision with only Spencer-owned gaps, sweep 45/45, E041-E048 arc
  settled on three sizes. Loop enters idle cadence; wrap-up digest in the
  final hours.

## Iteration 28 (2026-08-01 ~16:30) — wrap-up delivered early

- All work complete ~48 h before the hard stop; wrap-up digest written
  and delivered now (research/loop2_wrapup.md) rather than sitting idle.
  Loop stays on light heartbeats through the window to catch any new
  Spencer directives; hard stop unchanged (2026-08-03 ~19:12 JST).
