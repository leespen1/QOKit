# Research Next Steps

<!-- Priority: P0 = do next, P1 = do soon, P2 = backlog, DONE = completed -->
<!-- When completing an item, move it to DONE with a link to the log entry -->

## Research Focus — Updated 2026-04-02

**Paper thesis (revised):** The homogeneous parameter-setting heuristic, when
paired with sampling-based homodist estimation, produces competitive QAOA
parameters for Erdős-Rényi graphs across the full p_edge spectrum. Sparse ER
is the strongest testbed, where both PaperProxy (analytical) and SampN+EmpP
(empirical) beat Transfer. For non-ER graphs, parametric proxies (Triangle,
Normal, Gaussian) fail, and the proxy heuristic does not outperform Transfer.

**Key findings (consolidated from all experiments):**

1. **Proxy beats Transfer on ER graphs** — HEADLINE
   - At p=1: SampN+EmpP is best. Advantage grows with n (+2-3% at n=22).
   - At p=3: PaperProxy(correct p_edge) is best for sparse ER (+3.5-3.9% over Transfer at n=22).
   - At p=3: SampN+EmpP beats Transfer on ER(0.5) at large n (+3.6% at n=22).
   - At p=5: Transfer still wins (proxy landscape overfits).

2. **Sparse ER is the strongest testbed**
   - Random baseline: ER(0.2)=0.60, ER(0.3)=0.65 vs ER(0.5)=0.75.
   - Larger headroom means proxy improvements are more meaningful.
   - PaperProxy formula works for any p_edge — natural extension.

3. **Non-ER is not viable for proxy heuristic** — CLOSED DIRECTION
   - SampN+EmpP loses to Transfer on ALL non-ER (BA/WS), regardless of
     sparsity. Gap worsens with n.
   - GaussianProxy (7 params) and smoothed homodist both perform at/below random.
   - TriangleProxy and NormalProxy also failed (prior results).
   - The homogeneity assumption breaks down for structured graphs.

4. **N/P consistency is essential** — CONFIRMED
   - Mixing N from one source with P from another is catastrophic.
   - PaperProxy succeeds because N and P come from the same analytical model.
   - SampN+EmpP succeeds because both come from the same instances.

5. **Parametric proxy fitting is a dead end**
   - Homodist MSE ≠ QAOA performance. Good shape fit doesn't imply good parameters.
   - Error amplification across layers kills parametric proxies.
   - The only thing that works: use the actual homodist (empirical or analytical).

6. **Theoretical framework: compression vs transition error** (2026-04-03)
   - Proxy error decomposes into: (a) compression error (2^n → m+1 lossy projection),
     (b) transition error (wrong N within compressed space). Homodist MSE measures only (b).
   - Non-ER fails because (a) dominates: structured graphs have high within-cost-class
     amplitude variance that no choice of N can fix.
   - MSE ≠ performance because MSE is the wrong loss function: landscape gradient
     fidelity ≠ distribution fidelity. PaperProxy works despite errors because its N/P
     errors are structurally correlated and partially cancel.
   - See entry `2026-04-03_compression-vs-transition-error` for full analysis.

## P0 — Immediate

- [ ] **SVD of n(x;d,c): low-rank characterization of proxy validity**
  The homogeneous approximation IS a rank-(m+1) approximation of the n(x;d,c)
  matrix (rows=bitstrings, cols=(d,c) pairs). Compute SVD for small ER/BA/WS
  graphs and plot singular value spectra. Prediction: ER has fast decay
  (effective rank ≈ m+1), BA/WS have slow decay (high effective rank). This
  reframes proxy failure as a spectral property, connects to established
  low-rank approximation theory (Eckart-Young), and provides an a priori
  predictor of proxy quality. Could be a strong paper figure and theoretical
  contribution. Also check whether SVD basis vectors align with cost-class
  indicators (they should for ER, not for BA/WS).

- [ ] **Null hypothesis check: proxy vs random QAOA parameters**
  Are the proxy-optimized parameters actually better than random (γ,β) draws?
  If random parameters achieve similar approximation ratios, the proxy isn't
  helping — the "advantage over Transfer" could just be that any reasonable
  parameters work at low depth. Test: sample many random (γ,β) from the same
  range as proxy grid, evaluate real QAOA, compare distribution to proxy-optimal.
  Also compare to: random grid search (best of K random parameter sets). This
  is essential for credibility — reviewers will ask.

- [ ] **Paper figure design and generation**
  Based on all results, the paper needs these key figures:
  1. Random baseline context figure (all configs)
  2. ER(p_edge) comparison: proxy vs transfer with random baseline shown
  3. Sparse ER at large n: the headline result
  4. Non-ER negative results with analysis
  5. Sampling-based homodist: accuracy and computational scaling

- [ ] **Consolidate ER(0.5) + ER(0.2-0.3) results into unified table**
  Combine the headline results from `sampled-homodist-very-large-n`,
  `er-varying-pedge`, and `er-sparse-large-n` into one comprehensive
  comparison across the full ER parameter space.

## P1 — Soon

- [ ] **SampN+EmpP + Refinement vs Transfer + Refinement at large n**
  Script ready: `julia/paper_figures/sampn_refinement_large_n.jl`. Run on
  MSU HPC. Tests whether SampN+Refine extends the proxy advantage to higher
  depths where raw proxy optimization falters.

- [x] **ER(p_edge=0.2-0.3) at p=5 with Transfer comparison** → DONE
  See entry `2026-04-02_er-sparse-high-depth`. PaperProxy beats Transfer
  at p=5 on sparse ER(0.2) at n=20 (+0.8%). Proxy overfitting is ER(0.5)-specific.

- [x] **Sensitivity to number of homodist instances on sparse ER** → DONE
  See entry `2026-04-02_instance-count-sparse-er`. 5 instances sufficient at
  p=1, 10 recommended at p=3 for sparse ER. Not dramatically more than ER(0.5).

- [ ] **Literature review: low-rank methods in QAOA**
  The homogeneous approximation is effectively a low-rank projection of the full
  2^n state into an (m+1)-dimensional cost-class subspace. Search for prior work
  on: (a) explicit low-rank / tensor-network / matrix-product-state approaches to
  QAOA simulation or parameter setting, (b) compressed or reduced-basis methods
  for variational quantum algorithms, (c) any connection between spectral properties
  of problem Hamiltonians and QAOA performance. Compare to our SVD/spectral
  framing of the homogeneous heuristic — is this connection already known, or is
  reinterpreting the proxy as a low-rank approximation a novel theoretical lens?
  Key search terms: "low-rank QAOA", "tensor network QAOA", "QAOA compression",
  "reduced basis VQA", "QAOA symmetry reduction". Check arxiv 2023-2026.

- [ ] **Paper writing: abstract and introduction draft**
  With all experimental results in hand, draft the paper structure.

## P2 — Backlog

- [ ] **Measure within-cost-class amplitude variance (ER vs non-ER)**
  Testable prediction from compression error theory: compute Var[ψ(x) | c(x)=c']
  after each QAOA layer. Should be higher for BA/WS than ER, and should grow with
  depth. Would provide quantitative evidence for why the proxy fails on non-ER.
  This is a good paper figure if the effect is clear.

- [ ] **Per-instance homogeneity vs proxy performance correlation**
  Among ER instances, those with higher within-cost-class variance should have
  worse proxy performance. Would confirm compression error as the binding constraint.

- [ ] **ER with p_edge < 0.1 (very sparse, tree-like)**
  At very low p_edge, ER graphs become tree-like. Do proxy methods still
  work? This is theoretically interesting but may not be practical.

- [ ] **PaperProxy failure analysis at p_edge=0.7**
  Why does the analytical formula break at high density? Is it the binomial
  P(c') assumption or the multinomial N(c';d,c) formula? Understanding this
  could guide extensions to other graph families.

## DONE

- [x] **Replicate paper figures 2-7** (2026-03-17)
  Result: All 6 figures reproduced with small parameters. See `ClaudeResearchTasks/Task1/research_task1_final.md`

- [x] **Investigate proxy performance on BA/WS graphs** (2026-03-17)
  Result: PaperProxy w/ p_eff gets good Pearson but falls behind transfer at depth>1. Cross-type transfer works. BA/WS more homogeneous than ER. See entry `2026-03-17_smallworld-proxy-performance`

- [x] **Fit TriangleProxy to BA/WS empirical N(c';d,c)** (2026-04-01)
  Result: NEGATIVE. Better homodist fit but worse QAOA. N/P consistency matters more than fit quality. See entry `2026-04-01_fitted-triangle-proxy-evaluation`.

- [x] **Test fitted N(c';d,c) with consistent P(c')** (2026-04-01)
  Result: EmpN+EmpP best non-ER proxy. Fitted+BinomialP catastrophic. N/P consistency principle confirmed. See entry `2026-04-01_np-consistency-investigation`.

- [x] **Sweep BA/WS graph parameters** (2026-04-01)
  Result: PaperProxy fails at p_eff>0.5. Transfer robust across all configs. See entry `2026-04-01_parameter-sweep-robustness`.

- [x] **Scale to n=14-16** (2026-04-01)
  Result: Trends stable. Homodist feasible to n≈18 on CPU. See entry `2026-04-01_scaling-investigation`.

- [x] **Depth scaling p=1-8** (2026-04-01)
  Result: Transfer monotonic to 0.94+. Proxy oscillates. Gap ~0.08. See entry `2026-04-01_depth-scaling`.

- [x] **Hybrid proxy warmstart** (2026-04-01)
  Result: Warmstart+Refine beats Transfer by +0.01-0.03. See entry `2026-04-01_hybrid-warmstart`.

- [x] **Transfer+Refine vs Proxy+Refine** (2026-04-01)
  Result: Transfer+Refine ≥ Proxy+Refine. Proxy warmstart not uniquely valuable. See entry `2026-04-01_transfer-vs-proxy-refine`.

- [x] **Cross-type transfer at higher depths** (2026-04-01)
  Result: Works at all depths (penalty ≤ 0.016). BA/WS sources sometimes better than ER. See entry `2026-04-01_cross-type-transfer-depth`.

- [x] **Sampling-based homodist estimation** (2026-04-02)
  Result: S=10 samples/cost reproduces exact homodist. Speedup: 391x at n=18. See entry `2026-04-02_sampling-homodist-estimation`.

- [x] **Systematic TriangleProxy/NormalProxy evaluation** (2026-04-02)
  Result: NEGATIVE for alt proxies. EmpN+EmpP and SampN+EmpP tied as best non-analytical. See entry `2026-04-02_systematic-proxy-evaluation`.

- [x] **Why PaperProxy+p_eff outperforms empirical on non-ER** (2026-04-02)
  Result: Regularization effect. Low p_eff=0.1-0.2 best. See entry `2026-04-02_paperproxy-advantage-non-er`.

- [x] **SampN+EmpP on ER at n=14-18** (2026-04-02)
  Result: Beats Transfer on ER at n>=16. Degrades on non-ER. See entry `2026-04-02_sampled-homodist-large-n`.

- [x] **SampN+EmpP on ER at n=20-24** (2026-04-02)
  Result: HEADLINE. Beats Transfer by 3.6% at n=22,p=3. See entry `2026-04-02_sampled-homodist-very-large-n`.

- [x] **PaperProxy with low p_eff at larger n** (2026-04-02)
  Result: Inconsistent at high depth. Not reliable. See entry `2026-04-02_low-peff-large-n`.

- [x] **Dense non-ER investigation** (2026-04-02)
  Result: SampN+EmpP is only viable proxy for dense non-ER. See entry `2026-04-02_dense-nonER-investigation`.

- [x] **Higher depth (p=5) for SampN+EmpP at large n** (2026-04-02)
  Result: Transfer wins at p=5. Proxy landscape overfits. See entry `2026-04-02_high-depth-large-n-er`.

- [x] **Sample count sensitivity (S)** (2026-04-02)
  Result: S=5-10 sufficient. See entry `2026-04-02_sample-count-sensitivity`.

- [x] **Computational scaling analysis** (2026-04-02)
  Result: SampN+EmpP 67s at n=22. See entry `2026-04-02_computational-scaling`.

- [x] **Finer grid at p=5** (2026-04-02)
  Result: NEGATIVE. Finer grid worsens SampN+EmpP. Proxy overfits. See entry `2026-04-02_finer-grid-p5`.

- [x] **Multi-instance count sensitivity** (2026-04-02)
  Result: EmpN+EmpP saturates at ~10 instances. See `2026-04-02_paperproxy-advantage-non-er`.

- [x] **NormalProxy comparison** (2026-04-02)
  Result: Consistently worst. See `2026-04-02_systematic-proxy-evaluation`.

- [x] **Random baseline approximation ratios** (2026-04-02)
  Result: ER(0.5)=0.75, BA=0.62, WS=0.63, ER(0.1)=0.53. Biased partitions don't help. See entry `2026-04-02_random-baselines`.

- [x] **ER with p!=0.5 comparison** (2026-04-02)
  Result: PaperProxy strong for sparse ER(0.1-0.3), catastrophic at p=0.7. SampN+EmpP most robust at p=1. See entry `2026-04-02_er-varying-pedge`.

- [x] **Gaussian blob proxy** (2026-04-02)
  Result: NEGATIVE. Performs at/below random. Parametric fitting is dead end. See entry `2026-04-02_gaussian-proxy-investigation`.

- [x] **Smoothed homodist / spline proxy** (2026-04-02)
  Result: NEGATIVE. All smoothing destroys essential structure. See entry `2026-04-02_smoothed-homodist-investigation`.

- [x] **Sparse non-ER investigation** (2026-04-02)
  Result: SampN+EmpP loses to Transfer on ALL non-ER. Proxy is ER-specific. See entry `2026-04-02_sparse-nonER-investigation`.

- [x] **Sparse ER at large n (n=18-22)** (2026-04-02)
  Result: HEADLINE. Proxy beats Transfer on ER(0.2-0.3); PP best at p=3 (+3.5-3.9%); SampN best at p=1. See entry `2026-04-02_er-sparse-large-n`.

- [x] **Literature review: QAOA parameter-setting 2024-2026** (2026-04-02)
  See entry `2026-04-02_literature-review`. All contributions novel. 15 related papers, 10 must-cite.
