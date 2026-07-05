# Final referee verdict — qce2027_paper.tex (2026-07-05)

Compiled 9 pp draft, 7 figures, 3 tables; prior passes verified applied, not repeated.
E022–E024 checked against experiment READMEs; E023/E024 recomputed from raw CSVs.

## Verdict: weak accept

The exactness theorem plus leakage calculus stands, and E022–E024 supply what earlier
drafts lacked — a regime the proxy wins, located at a principled boundary, rechecked at
n=22–26, every traced number matching the records. Full accept is withheld only for
mechanical defects in the new material: a Limits paragraph the new experiments falsify,
boundary wording stronger than a 5-point one-family sweep supports, and missing N/SE for
5–10-instance cells judged under the paper's own 2xSE rule.

## Spot-checks of prior passes' worst items — PASS

- No "flat/stable in n" anywhere; l. 582 is the honest correction; Table II caption states
  the growth. All "recipe" mentions scoped ("Among proxy variants…", l. 713).

## New-material check (E022–E024)

Verified: E022 p=1 0.012–0.022 vs 0.031–0.099 (2–7x), p=3 "3 of 4 cells" with 0.021–0.140
vs 0.046–0.067; E023 α-grid; E024 0.015–0.025 vs 0.030–0.130 at α=1.5, reversal at α=3;
abstract/contribution-5/conclusion tell one story. Defects, by severity:
1. **Limits ll. 890–894 now false**: "claims are for unweighted MaxCut" and "only U[0,1]
   weights were tested" predate E021–E024 and contradict fig:alpha outright.
2. **"Locates the boundary precisely" (l. 793), bold "sits at the infinite-variance
   boundary" (ll. 798–799), caption "crossing at α=2"**: from the CSVs, 3-regular has
   already flipped to transfer at α=2 (0.0071±0.0023 vs 0.0336±0.0062, >4σ) — the sharp
   point is dense-ER-only. Soften to "consistent with the infinite-variance point α≈2,"
   report the 3-regular flip; 5 α values at one n, one weight family cannot "locate precisely."
3. **No N or SE in any weighted paragraph** (10/cell E022–E023; 5/cell E024) while §V-A
   promises 20–30 and a 2xSE rule. Two headline cells are ties under that rule: E024
   3-regular n=26 α=1.5 (0.0220±0.0115 vs 0.0299±0.0166) and E023 3-regular α=1.5
   best-of-3 (0.0119±0.0031 vs 0.0196±0.0059). Replace "ahead in every cell" with the
   defensible form (dense ER decisive, >4σ; 12/12 sign consistency); bands on fig:alpha.
4. Proxy "nearly α-independent (~0.02–0.03)": measured span 0.012–0.046, mild upward
   trend. Quote the real range.
5. Abstract drops the p=1 qualifier on "wins by 2–7x" (contribution 5 and conclusion keep it).
6. Contribution 5 omits the α≈2 boundary: the most quotable new result is in no contribution item.

## Three residual reviewer risks (one-sentence mitigations)

1. **Pareto-only**: infinite variance vs heavy-tails-per-se undiscriminated — add to
   Limits: "measured for Pareto tails; a finite-variance heavy-tailed law (lognormal)
   would separate the mechanisms and is untested."
2. **Beyond MaxCut**: the whole map is MaxCut while sud2024 covers k-XOR/SAT — add:
   "Thms 1–3 hold for any partition; only Cor. 4 and the map are MaxCut-measured; a
   Max-k-XOR replication is the immediate next experiment."
3. **Noiseless selection**: add "all methods select angles offline; hardware noise enters
   only through schedule robustness, common to every method compared."

## Length/figures (10-page limit with references)

9 pp + refs is tight. Claim-free cuts: fig:bound right panel (repeats fig:anatomy's
density ordering); fig:ranking → panel (a) plus one ρ-fade sentence; Lemmas 4–5 + Remark
→ 6-line remark (~0.3 pp; τ²/m detail is repo material); merge the duplicated
trajectory-PCA paragraph (§V-B end) into the Discussion E016 paragraph; fig:binned to one
column. Keep the glossary, fig:alpha, fig:directions, and Tables II–III.
