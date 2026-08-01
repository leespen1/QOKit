# Loop 2 wrap-up digest (2026-07-31 → 2026-08-01; window open to 08-03)

*All planned work finished ~48 hours early. One page; everything links to
the journal (loop_journal.md, iterations 1-26), an experiment, or a commit.*

## The one-paragraph version

The two diverged paper histories are reconciled (no history lost), the
manuscripts are submission-shaped, and the science advanced twice: the
noise-regularization arc (E041-E044: sampled N + normalized objective beats
even the exact compression, because structured sampling noise disrupts the
winner's curse) and the polynomial-time N estimator (E046-E048: the paper's
open question closes YES at low depth). Every number in both manuscripts
traces to a record, and every record regenerates. A fresh referee panel
answers your journal question: **yes — minor-revision equivalent at QINP
(Springer) tier; no further research is required for worthiness.**

## What happened, in order

1. **Reconciliation.** The paper repo had two diverged lines (the July
   editorial line was never pulled by loop 1). Resolved: July line
   canonical, loop-1 science ported in, histories joined, experiment
   numbering split (E017-E032 July / E033-E048 loop; map in STATUS).
   Two manuscripts now: `qinp_paper.tex` (journal flagship, 21 pp) and
   `qce2027_paper.tex` (conference cut, exactly 10 pp — your target).
2. **Verification gauntlet.** Dual audit (numbers vs CSVs + full proof
   re-derivation: no broken theorems, 21 precision fixes), a four-referee
   gauntlet (11 fixes incl. two rendering blockers), and a 45/45
   reproduction sweep (17 full reruns byte-identical, rest recompute from
   committed data, zero failures) — reproduction_sweep_2026-08-01.md.
3. **Science.** E041-E044: the composed recipe is the best proxy variant
   at p<=3, stable n=12-18; the mechanism is structured sampling noise
   (less data gives better parameters). E045 (negative): margin
   subtraction cannot fix the certificate's slack; selection-aware regret
   theory remains the sharply-delimited open problem. E046-E048: a
   Wang-Landau-based polynomial estimator matches the enumeration sampler
   at p=1 up to n=18, trails at depth (0.005/0.02/0.026 at n=14/16/18),
   estimates N in ~1.3 s at n=24 — in both manuscripts, prototype-grade
   scoping, with the "advantage grows with n" overclaim caught and
   corrected by E048 before it settled.
4. **Documents.** paper_explainer.md rebuilt against the flagship (10.3k
   words, every number provenance-tagged); journal_readiness.md refreshed;
   front matter done (abstract 197 words, keywords, declarations skeleton,
   all 35 references verified against publication records); all figures
   vector; final coherence read passed.

## Needs you (nothing else blocks submission)

1. **Declarations**: author list, funding, competing interests,
   contributions (placeholders marked TBD in qinp_paper.tex), and the
   frozen Zenodo/DOI archive.
2. **One framing choice**: the flagship now says the estimator is
   validated against truth to n=18 (matching its own numbers); E047's
   README frames n>=16 as "beyond validation." Pick a framing if you
   prefer the stricter one.
3. FYI: a 45-agent parallel workflow burned a session limit in one minute
   with zero output early on (your flag was right); everything after ran
   sequential/inline and the lesson is saved to memory.

## Where things live

STATUS.md (one page, current) · loop_journal.md Loop 2 (audit trail) ·
journal_readiness.md (verdict) · reproduction_sweep_2026-08-01.md ·
experiments 041-048 (each with Question/Answer/Method/Repro) · both
manuscripts on the paper repo's ClaudeResearch, all pushed.
