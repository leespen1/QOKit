# Research Log Index

<!-- One line per entry. Claude reads this to check novelty before logging new results. -->

| Date | File | Section | Summary |
|------|------|---------|---------|
| 2026-03-17 | smallworld-proxy-performance | parameter-transfer | PaperProxy w/ p_eff gets r>=0.98 Pearson on BA/WS but falls behind at depth>1; cross-type transfer works within 0.005 |
| 2026-04-01 | fitted-triangle-proxy-evaluation | approximation-ratios | FittedTriangleProxy achieves 10x lower homodist MSE but worse QAOA performance; internal N/P consistency matters more than fit quality |
| 2026-04-01 | parameter-sweep-robustness | approximation-ratios | PaperProxy catastrophically fails for dense non-ER graphs (p_eff>0.5); Transfer robust across all BA/WS parameter configurations |
| 2026-04-01 | np-consistency-investigation | methodology | EmpiricalN+EmpiricalP is best proxy for non-ER; fitted N + binomial P catastrophic; N/P consistency principle confirmed |
| 2026-04-01 | scaling-investigation | scalability | Trends stable n=12→16; Transfer best; EmpN+EmpP beats PaperProxy for ER; homodist feasible up to n≈18 on CPU |
| 2026-04-01 | depth-scaling | scalability | Transfer monotonic to 0.94+ at p=8; PaperProxy oscillates; gap plateaus ~0.08; all graph types behave similarly |
| 2026-04-01 | hybrid-warmstart | methodology | Proxy+Refine beats Transfer (no refine) by 0.01-0.03 at p≥3; warmstart valuable over random |
| 2026-04-01 | transfer-vs-proxy-refine | methodology | Transfer+Refine ≥ Proxy+Refine; proxy warmstart not uniquely valuable; refinement on target is the key improvement |
| 2026-04-01 | cross-type-transfer-depth | parameter-transfer | Cross-type transfer works at all depths (penalty ≤ 0.016); BA/WS sources sometimes better than ER even for ER targets |
| 2026-04-02 | sampling-homodist-estimation | methodology | Sampling-based homodist: S=10 samples/cost reproduces exact; 16-20x speedup at n=12; QAOA params robust to sampling noise |
| 2026-04-02 | systematic-proxy-evaluation | approximation-ratios | Tri+TriP degrades at p>1; NormalProxy worst; EmpN+EmpP≈SampN+EmpP; PaperProxy+p_eff surprisingly strong on non-ER |
| 2026-04-02 | paperproxy-advantage-non-er | discussion | PaperProxy advantage is regularization (low p_eff=0.1-0.2 best); more instances don't help EmpN+EmpP; P source irrelevant |
| 2026-04-02 | sampled-homodist-large-n | scalability | SampN+EmpP beats Transfer on ER at n>=16 (0.865 vs 0.846 at n=18,p=3); degrades on non-ER; computation practical (~1.5s at n=18) |
| 2026-04-02 | sampled-homodist-very-large-n | scalability | HEADLINE: SampN+EmpP beats Transfer by 3.6% on ER at n=22,p=3 (0.881 vs 0.844); advantage grows with n; ~1min compute at n=22 |
| 2026-04-02 | low-peff-large-n | scalability | PP(0.2) sometimes helps at p=3 on WS (+3.3%) but inconsistent at high depth; not a reliable recommendation |
| 2026-04-02 | dense-nonER-investigation | approximation-ratios | SampN+EmpP is only viable proxy for dense non-ER (p_eff>0.5); beats PaperProxy by 15% on BA(m=4); Transfer still best |
| 2026-04-02 | high-depth-large-n-er | scalability | SampN+EmpP beats Transfer at p=1,3 (n>=18) but Transfer wins at p=5 everywhere; advantage is depth-dependent |
| 2026-04-02 | sample-count-sensitivity | methodology | S=5-10 samples/cost is sufficient; S=5-100 produce identical QAOA performance; recommend S=10 as default |
| 2026-04-02 | computational-scaling | scalability | SampN+EmpP: 5.6s at n=18, 67s at n=22; refinement adds 34s at n=20; Transfer is free at deployment |
| 2026-04-02 | finer-grid-p5 | discussion | Finer proxy grid makes SampN+EmpP WORSE at p=5; proxy landscape overfits at high depth; motivates warmstart+refine |
| 2026-04-02 | literature-review | discussion | All 3 contributions novel; no prior sampling-based homodist; 15 related papers identified; 10 must-cite |
| 2026-04-02 | random-baselines | methodology | Random baseline: ER(0.5)=0.75, BA=0.62, WS=0.63, ER(0.1)=0.53; biased partitions don't help; sparse graphs have most headroom |
| 2026-04-02 | er-varying-pedge | approximation-ratios | PaperProxy strong for sparse ER(0.1-0.3); catastrophic at p=0.7; SampN+EmpP most robust at p=1; Transfer wins p=3 sparse |
| 2026-04-02 | gaussian-proxy-investigation | approximation-ratios | NEGATIVE: GaussianProxy (7 params) performs at/below random despite good homodist fit; parametric proxy fitting is a dead end |
| 2026-04-02 | smoothed-homodist-investigation | approximation-ratios | NEGATIVE: Smoothing homodist always worse than raw; even σ=0.5 destroys useful structure; fine detail is essential not noise |
| 2026-04-02 | sparse-nonER-investigation | approximation-ratios | SampN+EmpP loses to Transfer on ALL non-ER (sparse or dense); gap worsens with n; proxy is ER-specific |
| 2026-04-02 | er-sparse-large-n | approximation-ratios | HEADLINE: Proxy beats Transfer on sparse ER(0.2-0.3) at n=18-22; PP best at p=3 (+3.5-3.9%); SampN best at p=1 (+2-3%) |
| 2026-04-02 | er-sparse-high-depth | scalability | PaperProxy beats Transfer at p=5 on sparse ER(0.2) n=20 (+0.8%); proxy overfitting is ER(0.5)-specific not intrinsic |
| 2026-04-02 | instance-count-sparse-er | methodology | 5 instances sufficient at p=1; 10 recommended at p=3 for sparse ER; 20 provides minimal additional benefit |
| 2026-04-03 | compression-vs-transition-error | discussion | Proxy error = compression (2^n→m+1 lossy projection) + transition (wrong N). Non-ER fails because compression error dominates; MSE≠performance because MSE is wrong loss function |
