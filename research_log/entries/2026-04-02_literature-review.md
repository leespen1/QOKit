# Literature Review: QAOA Parameter-Setting Methods (2024-2026)

**Date:** 2026-04-02
**Paper section:** discussion
**Tags:** literature-review, novelty-check, related-work
**Status:** complete

## Motivation
Confirm novelty of our three contributions before IEEE Quantum Week submission.
Identify related work to cite and check whether sampling-based homodist has been
scooped. Prompted by P0 item in next_steps.md.

## Setup
- **Method:** Systematic web search of Google Scholar, arXiv, IEEE Xplore
- **Queries:** Citations of Sud et al. (2024), "QAOA parameter" 2024-2026,
  "QAOA sampling parameter", "QAOA parameter transfer", "QAOA Barabasi-Albert",
  "IEEE Quantum Week 2025 QAOA", "QAOA parameter landscape"

## Key Findings

### Novelty Assessment

| Contribution | Status | Closest Related Work |
|---|---|---|
| **Sampling-based homodist estimation** | **Suggested by Sud et al. but never implemented** | Sud et al. (2024) §Discussion explicitly suggest "estimate the distributions via Monte Carlo sampling of bitstrings and their costs for the given instance" as a future direction. No one has actually done it. Our contribution: implementation, validation (S=10 sufficient), scaling analysis, and the result that SampN+EmpP beats transfer at large n. |
| **N/P consistency principle** | **Novel** | Sud et al. use consistent N+P but don't identify it as a principle or show failure modes. |
| **PaperProxy as regularizer on non-ER** | **Novel** | Parameter transfer to non-ER exists but via optimization transfer, not proxy methods. |

### Direct Follow-ups to Sud et al. (Phys. Rev. Research 6, 023171, 2024)

1. **He, Shaydulin et al.** "Parameter Setting Heuristics Make QAOA Suitable for the Early Fault-Tolerant Era," ICCAD 2024 (arXiv:2408.09538). Survey of parameter-setting advances. Does NOT use homogeneous proxy. Complementary.
2. **Sureshbabu et al.** "Parameter Setting in QAOA of Weighted Problems," Quantum 8, 1231 (2024) (arXiv:2305.15201). Extends to weighted MaxCut. Analytical p=1 parameters. Does not use proxy for non-ER.
3. **Eichenseher et al.** "Pattern or Not? QAOA Parameter Heuristics and Potentials of Parsimony," arXiv:2510.08153 (2025). Parameters deviate from patterns at high depth. Does not use proxy.

### Parameter Transfer at Scale (Active Area)

4. **Galda, Liu et al.** "Transferability of optimal QAOA parameters between random graphs," IEEE QCE 2021 (arXiv:2106.07531). Foundational transfer work. 6→64 nodes with <1% loss.
5. **Montanez-Barrera, Michielsen** "Toward a linear-ramp QAOA protocol," npj Quantum Information (2025) (arXiv:2405.09169). LR-QAOA scales to 42 qubits/400 layers. Real QPU up to 109 qubits. Complementary — we use proxy to *optimize* ramp parameters.
6. **Dehn et al.** "Extrapolation method to optimize linear-ramp QAOA parameters," arXiv:2504.08577 (2025). Extrapolates LR-QAOA params from small→large. Alternative to our approach.
7. **Nguyen et al.** "Cross-Problem Parameter Transfer in QAOA: A Machine Learning Approach," arXiv:2504.10733 (2025). ML-based cross-problem transfer with graph embeddings.
8. **Hao et al.** "End-to-End Protocol for High-Quality QAOA Parameters with Few Shots," Phys. Rev. Research (2025) (arXiv:2408.00557). Transfer + fine-tuning on hardware. Up to 32 qubits on trapped-ion.
9. **arXiv:2601.15760 (2026).** Parameter transfer on 3-Regular, ER, and BA graphs. 98.88% optimal with 8x speedup.

### Non-ER Graph Studies

10. **Katial et al.** "On the Instance Dependence of Parameter Initialization for QAOA," INFORMS JoC 37(1) (2024) (arXiv:2401.08142). Instance space analysis across diverse graph types. Does not use proxy.
11. **Sang et al.** "Landscape-Similarity-Guided Optimization in QAOA," arXiv:2602.21689 (2026). Landscape universality across graph types. Proposes DO-QAOA.

### ML/Surrogate Approaches

12. **QAOA-GPT**, arXiv:2504.16350 (2025). Transformer generates QAOA circuits.
13. **QSeer**, arXiv:2505.06810 (2025). GNN-based parameter prediction.
14. **Deep-Circuit QAOA**, Quantum journal (2025) (arXiv:2210.12406). Classical performance indicator for deep QAOA.

### Comprehensive Reviews

15. **Blekos et al.** "A Review on Quantum Approximate Optimization Algorithm and its Variants," Physics Reports 1068 (2024) (arXiv:2306.09198). Main QAOA review paper. Must cite.

### IEEE Quantum Week 2025

No paper in QCE25 program addresses sampling-based homodist or homogeneous proxy for non-ER graphs. QAOA-GPT was presented there.

## Must-Cite Papers (Top 10)

1. Sud et al. (2024) — baseline paper
2. Blekos et al. (2024) — comprehensive QAOA review
3. He, Shaydulin et al. (2024) — parameter heuristics survey
4. Sureshbabu et al. (2024) — weighted MaxCut parameters
5. Galda, Liu et al. (2021) — parameter transferability
6. Montanez-Barrera, Michielsen (2025) — linear ramp scaling
7. Katial et al. (2024) — instance dependence
8. Hao et al. (2025) — end-to-end few-shot protocol
9. Eichenseher et al. (2025) — parameter parsimony
10. Sang et al. (2026) — landscape similarity

## Significance

**All three contributions are novel as of April 2026**, but with an important nuance:
Sud et al. explicitly suggest MC sampling of N(c';d,c) as a future direction in their
Discussion section. Our contribution is not the idea itself but its **realization and
validation**: we show S=10 samples/cost suffices, characterize computational scaling,
and demonstrate that SampN+EmpP beats transfer at large n on ER — a result that was
not obvious from the suggestion alone. The paper should frame this as "realizing and
extending the future direction proposed by Sud et al." rather than claiming the
sampling idea as wholly original. The N/P consistency principle and non-ER application
remain fully novel insights.

Note: The G-RIPS Sendai 2024 final report is publicly accessible and describes our
triangle/normal proxy work. The IEEE paper should reference it as prior work by us.

## Next Steps Arising
- [ ] Compile BibTeX entries for the 10 must-cite papers
- [ ] Position paper narrative: proxy heuristic for arbitrary graphs at scale
