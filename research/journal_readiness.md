# Journal readiness: adversarial verdict for the reconciled QINP manuscript

*Refreshed 2026-07-31 against `qinp_paper.tex` (single-column, ~20 pp,
Springer Quantum Information Processing target; paper repo branch
ClaudeResearch, head `1e72f46`). Method: three independent referee-style
reviews (novelty/significance, technical soundness, presentation/claims)
calibrated to QINP and comparable venues (EPJ Quantum Technology, Physica A,
IEEE TQE at the generous end). Supersedes the 2026-07-22 review, which was
written against the pre-reconciliation manuscript and an IEEE TQE/QST
target; that verdict is preserved in the history note at the bottom.*

## Direct answer to Spencer's question

**Yes: the manuscript is suitable, in both quality and novelty, for a
low-to-mid journal, and QINP specifically is a good fit.** No further
research questions need to be answered before the paper is worthy of
submission. What remains is administrative and editorial (declarations,
frozen data archive, figure format, one internal-consistency decision about
the E041-E043 arc), closable in days, not weeks.

The two candidate "more research first" items are optional:

- **E041-E043 (sampled-N noise regularization)**: not required. The paper is
  coherent and honest without it. But it is the single strongest available
  upgrade, and it creates a mild consistency obligation (see Gap 2): the
  paper's recommended recipe (sampled N with the unnormalized objective at
  p=1, normalize only at depth) is now known internally to be dominated by
  sampled N + normalized at both depths (E041: p=1 regret 0.021 vs 0.031;
  p=3 0.024 vs 0.073), with the mechanism resolved (E042: noise is the
  active ingredient, monotone in S, only under normalization; E043: the
  sampler's in-constraint-set noise structure is essential at depth). The
  quarantine condition E041 set ("do not put in the paper until E042
  resolves it") is now met. Folding it in would sharpen, not change, the
  paper's thesis (selection effects, not pointwise accuracy, govern
  parameter quality) and is recommended but optional.
- **Larger-n crossover demo (n=24-28 pipeline win)**: not required. The
  paper already scopes the timing claim honestly (Sec. timing: at simulable
  sizes the proxy buys no wall-clock advantage; the crossover is stated as
  extrapolation). A demonstrated win would strengthen the "beyond
  simulation" residue but the paper no longer stakes a speedup claim, so
  its absence is not a soundness gap.

## Referee A: novelty and significance. Verdict: MINOR REVISION

The core identification (the homogeneous proxy is exactly
QAOA-layer-then-project when driven by the instance's empirical N) is
simple to prove but genuinely new, and everything downstream of it earns
its place: the leakage calculus as a cheap a posteriori instrument, the
MaxCut first-order cancellation with explicit V2 functional, the density
law reduced to graph invariants (codegrees, triangles, 4-cycles) via the
polynomial-time conditioning bound, the two-point regret certificate with
its signed anatomy, and the argmax-transfer thesis with its controlled
direction experiments. The paper is candid that the telescoping bound is a
Dirac-Frenkel analogue and that Sud et al. already observed leading-order
homogeneity; novelty is claimed only for the object and the sharpening,
which is the right posture and survives scrutiny against the cited
literature (Sud 2024, Khairy-era surrogate fitting, Shaffer 2023,
fixed-angle/transfer work, Burgholzer 2025 exact bisimulations).

Two significance caveats keep this at "minor" rather than "strong accept":
(1) the headline question is answered largely negatively for the proxy as a
parameter setter (transfer is at least as good in every fair regime), so
the lasting contribution is the calculus and the diagnosis, not a method;
a referee could ask whether that carries a 20 pp paper. It does at this
tier, because the negative results are precise, controlled, and correct a
published open comparison (the Sud high-depth head-to-head), but it caps
the ceiling. (2) The theory is exact only for the empirical N, which is
unavailable exactly where the proxy would matter; the paper concedes this
openly and routes the beyond-simulation regime to the analytical model's
argmax robustness, which is the honest but modest reading. Neither caveat
is fixable by more experiments; they are the paper's actual shape.
"Proxy = exact projection + measured anatomy + practitioner map" is a
publishable core at QINP tier.

## Referee B: technical soundness. Verdict: MINOR REVISION

Spot checks of the theory chain and the statistics against the experiment
records found no errors and one loose end:

- Theorem 1 (compression), Theorem 3 (variance identity), Lemma 4
  (neighbor sum): elementary proofs check out; machine verification to
  1e-15/1e-16 recorded (E001, E008, unit tests).
- Corollary 5 (cubic law): the epistemic status of the linear-in-depth step
  is flagged in the text exactly as it should be (prediction validated
  empirically, ratio 1.505 vs 1.5, per-instance 1.54 +- 0.26 over 70
  index-paired runs with the ensemble-mismatch caveat stated; E011, E033).
- Proposition 7 (quadratic conditioning bound): the L2-projection argument
  is valid; moment identities machine-verified to 1e-9 on 140 instances,
  capture 91-100% (E035); weighted extension re-verified (E037).
- Proposition 8 (two-point certificate): machine-checked on 280
  instance-depth pairs; honest tightness reporting (median 10-11x) and the
  E038 signed decomposition match the records. **Loose end**: the
  unnormalized-objective addendum term is marked "audit-verified
  algebraically, not part of E036's machine check." Machine-check it or
  say why not; a careful referee will ask.
- Statistics: bootstrap CIs on every headline rho, family/size demeaning,
  within-cell and cluster controls, the pre-committed E004 gate stated as
  such, binomial p-values on sign claims, SEs on all family means, grid
  ceilings validated against continuous refinement (E033, E034). This is
  referee-proof in structure. Residual soft spots, all disclosed in the
  text: the n=16 displacement-regret correlation is convention-dependent
  (0.84/0.76 normalized vs 0.25/0.47 unnormalized), the weighted and
  refinement comparisons run on 5-10 instances per cell, and the ranking
  claim is scoped to p=1. These would draw questions, not a major-revision
  demand, because the paper already states them.

No claim in the manuscript contradicts its experiment record. The one
forward-looking soundness risk is Gap 2 below (the paper's recipe vs the
team's own E041-E043 knowledge), which is a completeness/honesty matter,
not an error.

## Referee C: presentation and claims. Verdict: MINOR REVISION (mechanical)

Scoping honesty is exemplary: negative results reported in full, the
E022-E026 heavy-tail correction narrated in the open, limits section
specific, every number traced to an experiment ID. Length (~20 pp
single-column) is normal for QINP. Internal consistency of numbers between
abstract, contributions, body, and records checks out on spot inspection.

Mechanical gaps that must close before submission:

- Springer front/back matter is absent: author list (placeholder),
  acknowledgments (placeholder), and the required declarations (funding,
  competing interests, data availability, code availability, author
  contributions).
- The reproducibility footnote points to a mutable git branch; QINP data
  policy and basic archival hygiene want a frozen DOI archive (Zenodo) of
  the experiment records and figure generators.
- All eight figures are included as PNG in `qinp_paper.tex`; the vector-PDF
  export pass applied earlier to the QCE cut needs to be pointed at this
  build, with fonts checked at single-column width.
- The abstract is a single dense ~230-word paragraph; fine for arXiv,
  worth one compression pass for the journal. Header hygiene: "Draft, not
  for distribution" date line, keywords/MSC codes missing, sn-jnl.cls swap
  is planned and fine to defer.

## Overall verdict

**MINOR REVISION equivalent: the science is submission-ready for QINP; the
package is not yet.** All three referees land at minor. Nothing on the list
requires new experiments or theory. Rank order of what must close versus
what would strengthen:

### Must close before submission (ranked)

1. **Declarations and archive.** Author list, acknowledgments, funding,
   data/code availability with a frozen Zenodo DOI replacing the mutable
   branch pointer. Blocking and Spencer-owned.
2. **The E041-E043 consistency decision.** Either fold the arc in (the
   E041 quarantine condition is met; it upgrades the practical section and
   the recipe changes from "normalize at depth" to "sampled N + normalized
   at both depths") or add a one-sentence forward pointer so the published
   recipe is not silently dominated by results already in the repo.
   Deciding is required; folding in is recommended.
3. **Figure production for the QINP build.** PNG to vector PDF, font and
   sizing check at single-column width.
4. **Certificate addendum parity.** Machine-check the
   unnormalized-objective extra term (or state plainly why the algebraic
   audit suffices), and give the n=16 convention-dependence one clarifying
   sentence where the 0.84/0.76 vs 0.25/0.47 split is reported.
5. **Front-matter pass.** Abstract compression, keywords/MSC, remove the
   draft watermark line, verify every bibliography entry against the
   published record (the loop has verified 36; re-verify arXiv-only items
   at submission time).

### Nice-to-have (not blocking, ranked by impact)

1. Fold in E041-E043 fully (if only the pointer is taken under Gap 2): the
   "less data gives better parameters" result is the most memorable single
   finding in the program and mechanistically complete.
2. A demonstrated pipeline win at n=24-28 (sampled N + proxy sweep vs GPU
   grid) to convert the extrapolated crossover into a measurement. Needs
   HPC; Spencer-owned.
3. Close the degree-3 residue of the conditioning correction, turning the
   density law fully into a theorem.
4. A polynomial-time Monte Carlo estimator of N (stated as open in the
   paper; a genuine follow-up paper, not a revision item).

## History

Prior verdict (2026-07-22, against the pre-reconciliation
`qce2027_paper.tex` at 12 pp, IEEE TQE/QST tier): morning review said major
revision on three decisive objections (significance ceiling from the
empirical-N wall; theory bounding state error while the question is argmax
displacement; missing uncertainty quantification and external baseline).
Same-day work closed them: E033 (statistics hardening; argmax claim
survives all controls, fidelity claim corrected), E034 (ceilings validated;
transfer baseline added and it beats the proxy), E035 (conditioning bound),
E036 (regret certificate + normalization rule). Evening re-review upgraded
to minor revision with only Spencer-owned placeholders blocking. All eight
ranked gaps from that review are closed except the HPC crossover demo
(item 4 there, nice-to-have item 2 here). The 2026-07-31 reconciliation
ported all of that into the QINP manuscript reviewed above.

## Gap status update (2026-07-31, loop 2 iteration 8)

- Gap 2 (E041-E043 consistency): CLOSED. The composition paragraph is in
  qinp_paper.tex section 6.1, scoped as best-proxy-variant-at-p<=3; the
  exact-vs-sampled normalization convention split is stated explicitly.
- Gap 3 (vector figures): CLOSED. All six QINP figures regenerated as
  vector PDFs from their committed generators (paper commit c2857c3).
- Gap 4 (addendum machine check + n=16 sentence): CLOSED. E036
  verify_addendum.jl asserts at 714 points; the n=16 convention split is
  stated in the manuscript's argmax paragraph.
- Remaining: gap 1 (declarations + frozen archive; Spencer-owned) and
  gap 5 (front-matter pass: abstract compression, keywords, watermark,
  bibliography re-verification).
