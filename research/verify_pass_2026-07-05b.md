# Verify pass 2026-07-05b — qce2027_paper.tex, E025–E031 passages + V₂ Lemma/Prop + Conclusion

Scope: passages tagged E025–E031, Lemma "Quartic split" / Prop "Conditioning correction",
Conclusion. Numbers recomputed from results_task*.csv (E025/026/028/029/030/031) and
checked against research/v2_conditioning_derivation.md.

## Findings

1. **Contradiction** (tex 827–829): "winning cleanly for $\alpha>2$" — E026 (eight lines
   later, tex 839–843) shows the proxy *winning* dense ER at α=2.2 (0.024 vs sources
   0.043–0.068). Fix: "winning cleanly at $\alpha\ge3$" (E023's grid had no point in (2,3)).
2. **Number-error** (tex 871–873): "the proxy still beats every \emph{single} weighted
   source ($0.05$--$0.11$ vs.\ $0.07$--$0.29$)" — in the lognormal ER p=10 cell the best
   source (0.0729) marginally undercuts both binned proxies (0.0731/0.0733): a dead tie,
   not a win. Also max source regret is 0.298 → "0.29" under-rounds. Fix: "beats or ties
   every single source ($0.05$--$0.11$ vs.\ $0.07$--$0.30$)".
3. **Number-error, minor** (tex 881–882): "wherever ${\sim}2000$ true evaluations are
   affordable" — E029's budget is ≤2000 *per stage*, two stages (ramp4 then full-2p);
   the 0.000–0.006 residual needs both, so up to ~4000. Fix: "${\sim}2000$ per
   refinement stage" or "a few thousand".
4. **Number-error, minor** (tex 887): "$0.80\to0.92$ lognormal" — measured 0.799→0.915
   (E030); 0.92 rounds up. Fix: "$0.80\to0.92$" → "$0.80\to0.91$" (Pareto
   $0.84\to0.88$ = 0.837→0.884 is fine; gain range +0.046…+0.116 → "+0.05–0.12" fine).
5. **Number-error** (tex 1004–1005): degree-≥3 residue "($4$--$6\%$ of it empirically)" —
   derivation doc: quadratic projection captures 94–99.7% per instance (family-mean
   residues 1.6–4.0%, ER(0.25) at 1.6%). Fix: "$1$--$6\%$" or "${\approx}4\%$".
6. **Contradiction, wording** (tex 67–70 abstract; 1023–1024 Conclusion): both keep the
   pre-E025 framing — the proxy's one won regime is "heavy-tailed weights" — while
   contribution 5 (tex 160–163) and the body (tex 838–844) place the boundary at weight
   *dispersion* (finite-variance lognormal σ≥1.5 is also proxy-won, E025, regret down to
   0.004). Fix: "high-dispersion (e.g.\ heavy-tailed) weights" in both places.
7. **Contradiction, minor** (tex 160–163 vs 843–844): contribution 5's dispersion clause
   omits the family-dependent threshold constant the body and E026 emphasize; and Limits
   (tex 968–974) inventories the crossover evidence as "a five-point $\alpha$ sweep …
   plus a two-$\alpha$ check … plus a lognormal $\sigma$ sweep" but omits the
   pre-registered three-point follow-up (E026, α∈{2.2,2.4,2.6}, 10/cell) that the
   headline "predictive" claim rests on. Add it to both.
8. **Scoping** (tex 876–899, 965–974): E029 ("grid artifact", residual 0.000–0.006) and
   E030 (polish +0.05–0.12; full-2p overfit 0.83) rest on 5 instances × 2 cells; E031's
   hard-corner numbers (0.77–0.80→0.90, best-of-two 0.92–0.95) on one 5-instance cell.
   Neither the body nor Limits states these counts (Limits gives counts for the α/σ
   sweeps only). Add "the depth-refinement and recipe checks use 5 instances per cell".

## Verified correct (no action)

- E025: transfer 0.16–0.35 at σ=2.5; proxy →0.004; flip σ≈1–1.5; CV²≈2–8; 10/cell.
- E026: dense-ER flip at α≈2.4; 3-regular below 2; fig:alpha caption numbers; 10/cell.
- E028: proxy 0.05–0.11 lognormal depth; Pareto "toward parity" properly grid-qualified;
  10/cell (not 5); no surviving unqualified "hard for every method" phrasing.
- E029: residual 0.000–0.006; 8⁴ ceilings 0.020–0.032 low (tex "0.02–0.03" OK).
- E031: 0.897→0.837 unweighted (0.8965→0.8367), single occurrence, guard story coherent;
  Pareto p=3 +0.05 (0.764/0.770→0.819); lognormal +0.093/+0.124; best-of-two
  0.919/0.951 → "0.92–0.95"; counts 10/10/5 per cell.
- Lemma quarticsplit + Prop condcorr: every identity, the V_q formula, κ(p_e), 94–99%,
  36τ²/m, Cov(T,S³) form all match v2_conditioning_derivation.md exactly.
