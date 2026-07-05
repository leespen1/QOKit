# Coherence pass: qce2027_paper.tex (referee, coherence-only)

Line numbers = current ClaudeResearch checkout of papers/OverleafPaper/qce2027_paper.tex.
Tags: [C] contradiction, [S] stale phrase, [R] ref error, [N] number mismatch, [G] contributions gap.

1. [C] L551, Table II caption: "Regret is flat in $n$" — contradicts the corrected body
   (L566-571: "grows only slowly with $n$ ... We previously read this as flat; the per-family
   table corrects that") and the table's own numbers (3-reg 0.037→0.055; every family monotone).
   Fix: caption → "Regret grows slowly with $n$ and roughly doubles from $p=1$ to $p=3$ ramps."
2. [C] L65 "is the best \emph{proxy} variant" (also L145, L701): unscoped, but E019 (L737-738)
   has exact compression (0.07–0.14) beating sampled (0.11–0.21) at p=10/20 — and that sentence
   itself counts exact compression among "every proxy variant". Fix: scope all three to p<=3.
3. [C] L817: "hardware parameter setting where every objective evaluation costs shots" listed as
   a proxy regime — contradicts L103-105 ("\emph{all} offline methods ... avoid that cost
   equally") and L725-727 ("Hardware does not change this"). Fix: cut the hardware item from the
   timing regime list (keep beyond-simulation sizes and GPU-less environments).
4. [C] L875 "the proxy wins none of the regimes we could test": the weighted regime (L754-766)
   was tested with no transfer baseline — research/experiments/021_weighted-transfer is
   "Pending" and uncited; abstract L66 "beats every proxy variant tested" likewise silently
   excludes the binned proxy. Fix: add a scope sentence in the weighted paragraph or Limits
   ("transfer not yet calibrated on weighted instances"); conclusion → "unweighted regimes".
5. [S] L97-98 "at every size we tested ($n\le20$)" (+ "by extrapolation ... $n\approx22$--$24$",
   L99-100): superseded by E019's GPU sweeps at n=22,26 (L741) and Limits' "parameter-quality
   comparisons to $n=26$" (L833). Fix: rescope to the timing runs; the crossover is now measured.
6. [S] L617 subsection title "Model error: robustness, one sharp failure, and a recipe": the
   recipe moved to Sec. sec:practical (L690-710). Fix: drop "and a recipe" from the title.
7. [R] L106 "(Sec.~\ref{sec:modelerror})" attached to "an empirical question we answer per
   regime": the per-regime map (transfer/sampled/analytical) is in sec:practical. Fix ref.
8. [R] L564 "the sampled-$N$ variant of Sec.~\ref{sec:modelerror}": sampled-N is defined in
   sec:practical (L692-698). Fix: ref sec:practical.
9. [R] L809 "the recipe of Sec.~\ref{sec:modelerror} is also the cheap option": the recipe lives
   in sec:practical (L701). Fix: ref sec:practical.
10. [R] L752 "the self-consistency effect of Sec.~\ref{sec:modelerror} again": that subsection
    never names or states a self-consistency effect (closest: L639-641, L686-688). Fix: name the
    effect there at first occurrence, or rephrase L752 without the back-reference.
11. [R] L700 "$\eta_F$ to median $3.2\%$": symbol never defined in the paper (E008-internal name
    for aggregate leakage). Fix: define at first use or write "aggregate leakage" in words.
12. [N] Bin count three ways: abstract L63-64 "${\sim}100$ bins", body L762 "$K\approx64$--$128$",
    Fig. 6 caption L775-776 "saturation by $K\approx64$". Fix: align; if p=1 saturates by 64 and
    p=3 by 128, say exactly that in the caption.
13. [N] Standard errors three ways: L474-475 "0.002–0.006", Table II caption L550-551
    "~0.001–0.006", L570 "the ~0.001 per-cell standard errors" — the growth claim leans on the
    smallest, and +0.005 fails the L476 "twice these" rule against 0.006. Fix: state regret-cell
    SE (~0.001, per generated table header) separately from AR-mean SE (0.002–0.006); use the
    right one at L570.
14. [N] L562 "0.031–0.055 at $n=16$--$18$" vs Figures/generated/regret_table.tex, whose minimum
    over n=16–18 at p=1 is 0.032 (ER(0.5), n=16). Fix: 0.031 → 0.032 (or requote from source).
15. [G] Contributions L108-151: the weighted-MaxCut/quantile-binning extension (abstract L61-64;
    full paragraph + two-column Fig. 6; E020) appears in no contribution item. Fix: add it —
    any-partition validity of Thms 2-3, lambda^2 ≈ lambda_struct^2 + O(1/K), O(K^2 n) binned
    proxy — to item 5 or as its own item.
16. [N] Symbol collision: Table III caption L785 "$K$ is the number of candidate schedules swept"
    vs K = number of bins throughout L758-777 and abstract L63. (The table's O(K m^2 p) also
    drops the n factor of the intro's "O(nm^2)-per-layer", L81.) Fix: rename the schedule count
    in Table III (e.g. J or "schedules").

Checked and found consistent (no action): "decouples at p=3" vs E019 (decoupling persists at
p=10/20, L737-741); sampled-N 0.011/0.016 (L143-144 vs L694-696); transfer dense-ER trend 0.04
at n=18 (L717-718) → 0.047→0.061 at n=22→26 (L743); slack 3–6x / worst ≤5.8x / median ≈4x;
exponents 0.98/1.92 and depth ratio 1.505 vs 1.5; weighted regret cross-checks 0.032 and 0.105
match Table II; all % E-tags E001–E020 match research/experiments/ directory contents; all \ref
labels resolve (items 7-9 are wrong-target, not dangling).
