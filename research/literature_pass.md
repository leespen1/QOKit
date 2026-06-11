# Literature pass on the claimed-new items (2026-06-11)

*Systematic check of the four claimed-new items in `theory_notes.md` against
the current literature (WebSearch + Semantic Scholar citation list of
arXiv:2211.09270, ~22 citing papers as of today). Verdicts: **still-novel**,
**needs-citation** (related, cite it, no preemption), or **overlaps** (none
found). Bib entries added to `papers/OverleafPaper/qce2027_paper.tex`.*

## Verdict summary

| Claim | Verdict |
|---|---|
| (a) proxy = exact orthogonal compression (Thm 1) | **still-novel** |
| (b) per-layer leakage + telescoping bound for QAOA (Thm 2) | **still-novel as applied**; technique standard (we already say so) |
| (c) f_d(β)-weighted within-class variance identity (Thm 3) | **still-novel** |
| (d) empirical program (leakage⇔fidelity⇔regret, density-not-ER-ness, error separation, argmax robustness/transfer) | **still-novel**; several needs-citation neighbors |

No citer of Sud et al. reinterprets the proxy as a projection/compression,
none gives error bounds for it, none fits parameterized shapes to N(c';d,c),
and none studies it on non-ER families. The citation list is dominated by
parameter-setting empirics, surveys, and applications.

## (a) Proxy as exact compression — checked against

- **Shaydulin et al. 2012.04713; Tsvelikhovskiy et al. 2309.13787** (already
  cited): *exact* symmetry subspaces. Ours is the approximate-invariance
  generalization; unchanged.
- **Allcock et al., "Analyzing QAOA: ansätze, symmetries, and Lie algebras"
  ([2410.05187](https://arxiv.org/abs/2410.05187))**: invariant subspaces of
  the generated Lie algebra for MaxCut — exact algebraic reduction, not
  cost-class lumping, no error quantification. Needs-citation alongside the
  symmetry line (optional).
- **"Hamiltonian-Guided Leverage Embedding" ([2606.07814](https://arxiv.org/abs/2606.07814))**:
  despite the "subspace compression for QAOA parameter estimation" title, it
  compresses *classical measurement-sample feature matrices* by leverage-score
  row sampling. Different object entirely; no cost classes, no dynamics
  projection. No threat; citation optional.
- **Quantum coarse-graining generalities ([1801.09770](https://arxiv.org/abs/1801.09770),
  Dávalos et al. PRA 2025)**: generic CG of quantum dynamics, no QAOA, no
  cost-class partition.

## (b) Telescoping leakage bound — checked against

- Projection/MOR machinery applied to quantum dynamics exists and is active:
  **"Model Order Reduction for Quantum Molecular Dynamics"
  ([2509.07340](https://arxiv.org/abs/2509.07340))**, reduced-order TDSE
  modeling with greedy a-posteriori estimators (Owolabi, Adv. Phys. Res.
  2026), measurement-adapted time-coarse-graining for open systems (PRX
  2025). None touches QAOA or partition (lumping) subspaces of {0,1}^n. Our
  framing — λ_ℓ measurable at one statevector layer + O(2^n) projection,
  empirically ~4× from tight — stands. Cite the MOR-for-quantum line as
  context if space allows.
- **Terminology collision to avoid**: "Mechanism of Efficacy in QAOA for
  Random k-SAT" ([2605.20288](https://arxiv.org/abs/2605.20288), 2026) uses
  "adiabatic leakage" = excitation out of the instantaneous adiabatic
  manifold. Unrelated to our subspace leakage; consider one disambiguating
  clause where we first define leakage.
- Symmetry-verification error mitigation ([2106.04410](https://arxiv.org/abs/2106.04410),
  [2204.05852](https://arxiv.org/abs/2204.05852)) projects *hardware-noisy
  states* onto symmetry sectors — different purpose; optional citation.

## (c) Weighted variance identity — nothing found

Searched Sud citers and QAOA-variance phrasings; nothing computes
class-size-weighted within-class variance of the f_d(β)-weighted combination
or links it to fit-norm choice. Still novel.

## (d) Empirical program — needs-citation neighbors (no preemption)

- **Krüger & Mauerer, "Out of the Loop: Structural Approximation of
  Optimisation Landscapes and non-Iterative Quantum Optimisation"
  ([2408.06493](https://arxiv.org/abs/2408.06493), Quantum 9, 1903 (2025))** —
  the closest paper in spirit to our argmax-transfer findings (E2.4/exp 012):
  approximates the p=1 QAOA landscape from solution-space structure,
  instance-independent but problem-specific, and proves parameter
  clustering. Complements rather than preempts: they build a good landscape
  surrogate; we diagnose *which error norms* of a surrogate predict
  parameter-setting regret (answer: argmax displacement, not entrywise MSE or
  amplitude error). **Cite in Related Work and §6.**
- **Sureshbabu et al., "Parameter Setting in Quantum Approximate Optimization
  of Weighted Problems" (Quantum 8, 1231 (2024))** — analytical rescaling
  rules for weighted instances. Cite with the transfer line.
- **Misra-Spieldenner et al., "Mean-Field Approximate Optimization Algorithm"
  ([2303.00329](https://arxiv.org/abs/2303.00329), PRX Quantum 4, 030335)** —
  classical spin-dynamics surrogate of QAOA; the other major
  classical-surrogate family. Cite in Related Work.
- **"Setting angles in QAOA at utility-scale" ([2606.05311](https://arxiv.org/abs/2606.05311),
  2026)** — fresh benchmark of angle-setting methods incl. small-to-large
  transfer at 100+ qubits; useful intro citation that parameter setting
  remains the practical bottleneck.
- **"Going off Pattern?" ([2510.08153](https://arxiv.org/abs/2510.08153))** —
  empirical: good parameters often deviate from expected patterns; supportive
  context, optional.
- Density-not-ER-ness: no published analogue found for compression/leakage
  vs graph density. Still ours.

## Search coverage caveats

- Semantic Scholar's citation list (~22 papers) and several WebSearch sweeps
  (projection/lumping/MOR/symmetry/surrogate phrasings, June 2026). Not a
  guarantee of exhaustiveness — very recent preprints and non-indexed venues
  can be missed; recommend one re-check at submission time (April 2027).
- Verdicts for (a)–(c) rest on abstracts; none was close enough to warrant a
  full-text audit.
