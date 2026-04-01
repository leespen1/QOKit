# Research Next Steps

<!-- Priority: P0 = do next, P1 = do soon, P2 = backlog, DONE = completed -->
<!-- When completing an item, move it to DONE with a link to the log entry -->

## P0 — Immediate

- [ ] **Fit TriangleProxy to BA/WS empirical N(c';d,c) and evaluate QAOA performance**
  Context: Smallworld investigation showed PaperProxy w/ p_eff falls behind transfer at depth>1 (entry 2026-03-17). TriangleProxy/NormalProxy fitted to non-ER data could do better since they can capture structure-dependent corrections.
  Approach: Generate many BA/WS instances (n=10, 50+ instances), compute empirical homodist, fit TriangleProxy via `sendai_opt.fit_proxy_to_real()`, then run figure6-style approx ratio comparison.
  Paper section: approximation-ratios

- [ ] **Sweep BA/WS graph parameters to test robustness**
  Context: Smallworld investigation tested only m_attach=2 (BA) and k=4,p_r=0.3 (WS). Need to know if results hold across parameter space.
  Approach: Vary m_attach in {1,2,3,4} for BA; vary (k,p_r) in {(2,0.1),(4,0.3),(6,0.5)} for WS. Run transfer + proxy comparison at p=1,2,3.
  Paper section: approximation-ratios

- [ ] **Scale to larger graphs (n=14-18) for non-ER types**
  Context: All current non-ER results are n<=12. Need to verify trends hold at paper-relevant sizes. GPU acceleration available for cost distribution computation.
  Approach: Run figure7-style depth scaling on BA/WS at n=14,16,18 using GPU-accelerated homodist.
  Paper section: scalability

## P1 — Soon

- [ ] **Depth scaling (p=1 to p=10) with linear ramp for all graph types**
  Context: Direction 3 — at higher depths the proxy landscape becomes more structured. Need to characterize the proxy-vs-transfer gap as a function of depth.
  Approach: Run figure7 script with BA/WS graph types, depths 1-10. Compare proxy and transfer approximation ratios.
  Paper section: scalability

- [ ] **Hybrid proxy warmstart + local refinement**
  Context: Direction 4 — use proxy to find good starting point, then refine with real QAOA + COBYLA. Could close the proxy-vs-transfer gap at higher depths.
  Approach: Run proxy optimization, take best params, run COBYLA on real QAOA starting from those params. Compare against pure transfer and pure proxy.
  Paper section: methodology

- [ ] **Statistical significance analysis**
  Context: Smallworld investigation noted lack of significance tests. Need proper error bars and p-values for the paper.
  Approach: Bootstrap confidence intervals or paired t-tests for approximation ratio comparisons. Increase instance count to 50-100.
  Paper section: methodology

## P2 — Backlog

- [ ] **Multi-instance averaged homodist for non-ER graphs**
  Context: Direction 5 — compute empirical N(c';d,c) averaged over many instances to get class-level smoothing without analytical formulas. BA/WS are more homogeneous (CV 2-3x lower), so this could work especially well.
  Approach: Average homodist over 100+ instances per graph type, use as proxy input.
  Paper section: methodology

- [ ] **Better fitted proxy shapes (TriangleProxy on ER)**
  Context: Direction 1 — fit TriangleProxy/NormalProxy to ER data averaged over many instances. Benchmark against PaperProxy to establish baseline before extending to non-ER.
  Paper section: distribution-analysis

- [ ] **NormalProxy comparison alongside TriangleProxy**
  Context: Two GRIPS proxy types exist but haven't been systematically compared on non-ER graphs.
  Approach: Run same experiments as TriangleProxy P0 item but with NormalProxy. Compare fit quality and QAOA performance.
  Paper section: approximation-ratios

- [ ] **Cross-type transfer at higher depths**
  Context: Smallworld investigation showed cross-type transfer works at p=1. Does this hold at p=2,3,5?
  Approach: Extend cross-type transfer experiment to higher depths using linear ramp grid.
  Paper section: parameter-transfer

## DONE

- [x] **Replicate paper figures 2-7** (2026-03-17)
  Result: All 6 figures reproduced with small parameters. See `ClaudeResearchTasks/Task1/research_task1_final.md`

- [x] **Investigate proxy performance on BA/WS graphs** (2026-03-17)
  Result: PaperProxy w/ p_eff gets good Pearson but falls behind transfer at depth>1. Cross-type transfer works. BA/WS more homogeneous than ER. See entry `2026-03-17_smallworld-proxy-performance`
