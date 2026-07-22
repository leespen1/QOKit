# Journal readiness: adversarial verdict and gap list

*Maintained by the improvement loop; this version 2026-07-22. Method: three
independent referee-style reviews of `qce2027_paper.tex` (post fold-in,
post audit, 12 pp), one per lens: novelty/significance, statistics/
methodology, claims/presentation. Full reports in the loop transcript;
everything actionable is distilled here.*

## Verdict

**Major revision at a low-to-mid journal (IEEE TQE, QST); not yet
submission-ready, but no referee expected the conclusions to reverse.**
Consensus strengths: a genuinely new and productive identification
(proxy = exact compression) of a well-cited heuristic; an unusually
rigorous, reproducible empirical program (fixed seeds, SEs, claim-to-
experiment mapping); honest negative results, one of which (the fitted-shape
paradox + argmax-transfer finding) is the strongest empirical content.

The three decisive objections:

1. **Significance ceiling (novelty referee).** The theory is exact only for
   the empirical N, which is unavailable in the regime where the proxy
   matters; the timing table concedes brute force matches the pipeline at
   every size tested, and the n≈22-24 crossover is extrapolated from one
   instance per size. The two findings that survive to the scalable regime
   (analytical-model argmax robustness; dense-ER artifact + norm
   certificate) should carry more of the framing.
2. **The theory bounds the wrong error (novelty referee, sharpest point).**
   The paper's own E010/E014 show leakage/fidelity decouple from regret
   beyond p=1; the apparatus predicts state error while the title question
   is answered by argmax displacement, about which the theory says nothing.
3. **Headline statistics lack uncertainty quantification (methodology
   referee).** No CIs on any of ~10 Spearman rhos; the pooled 140-instance
   correlations conflate between-family and within-family variation
   (Simpson's-type confound); the p=3 grid ceiling (8^4 linear ramps) is
   never validated against continuous optimization; no external baseline
   (published fixed angles / parameter transfer).

## What would flip the verdict (ranked by impact / cost)

| # | Gap | Fix | Cost | Status |
|---|---|---|---|---|
| 1 | Pooled rho confound | Family-demeaned (within-family) re-analysis of E012/E014 correlations + bootstrap CIs on every reported rho; state E004's pre-committed gate | Re-analysis of existing CSVs | **planned: E017** |
| 2 | Grid ceiling unvalidated | Continuous refinement (Nelder-Mead/BFGS) from grid argmax on a subsample; report ceiling gap; add a published-fixed-angle baseline row | CPU statevector at n<=14, minutes | **planned: E018** |
| 3 | Theory-regret gap | Any bound linking landscape sup-error (via leakage) to regret through peak flatness, even p=1-only | Theory, hard; highest impact | open |
| 4 | Crossover extrapolated | Demonstrate the pipeline actually winning at n=24-28 (sampled N + proxy sweep vs. GPU grid search) | Needs HPC GPU (local has none) | parked for Spencer/HPC |
| 5 | Density law empirical | Exact Var(E[T|c]) per ensemble (E015 reduced it to this) | Theory, "open and looks tractable" | open |
| 6 | Naked point estimates | 1.505 +- spread (E011 CSV); binomial p for 26/28 sign claim; single-instance hedge on timing caption | One evening | **planned: E017** |
| 7 | Fitted-shape target is internal practice | Either cite a published instance of entrywise N-fitting or present Table 3 explicitly as a cautionary internal replication | Literature check, small | open |
| 8 | Figure production | Re-export all six figures as vector PDFs, publication fonts, drop the "gate criterion" annotation, unify Gaussian/Normal naming | Re-run make_figure.jl scripts | open |

## Mechanical claim-discipline fixes (applied by the loop, 2026-07-22)

- "never vacuous" scoped (Fig. 2's extreme-ramp points exceed the trivial
  bound of 2).
- Abstract's 1.505 ratio now carries the small-ramps qualifier.
- "identifies the correct norm" softened to the state-error scope (the
  paper's own Sec. 5 dethrones every norm for parameter quality).
- Verification precision unified; abstract O(2^n) corrected to O(n 2^n);
  "trace to" softened to "connect to" for the conditioning cancellation.
- 9-vs-6 model-count mismatch between Fig. 5 caption and Table 3 reconciled.
- Note-to-self about re-running the literature pass removed from Limits.
- Self-description tics removed ("honest(ly)" x2, duplicated PCA
  disclaimer); E004 gate pre-commitment now stated in Sec. 5.3.
- Timing caption marks the crossover as a single-instance extrapolation.

## Submission-blocking placeholders (Spencer)

Author list; acknowledgments; frozen public archive (Zenodo/DOI) instead of
the mutable ClaudeResearch branch pointer; venue decision (file header says
IEEE TQE, filename says QCE 2027; QCE two-column will not fit Table 1 as
designed); verify the utilityscale2026 arXiv identifier before submission.

## Bottom line

The paper is one solid analysis pass (E017), one cheap validation experiment
(E018), and a figure-production pass away from a defensible major-revision
submission; the highest-impact remaining science (items 3-5) is what would
lift it from "solid contribution" toward "strong accept."
