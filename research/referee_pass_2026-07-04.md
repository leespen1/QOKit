# Referee pass — qce2027_paper.tex (2026-07-04)

*Unsympathetic-but-fair IEEE QCE referee, reviewing the current 8 pp IEEEtran
draft (post-E014/E017/E018 additions). Does not repeat the 2026-06-13 critique;
all its actioned items were verified as fixed. Line numbers refer to
`papers/OverleafPaper/qce2027_paper.tex` at paper-repo head.*

## Verdict

Accept-trajectory: the exactness theorem plus the leakage/variance error
calculus (Thms 1+3 — a new object under a known Dirac–Frenkel-type bound, with
an honest empirical anatomy) is a clear, well-foregrounded best case that sits
above the QCE bar, and the negative results are unusually credible. The one
thing that most needs fixing: the scale-stability headline — bolded "does not
grow with $n$" (l. 570), echoed in the abstract and contribution 4 — is
contradicted by the paper's own Table II, and the sampled-$N$ "recommended
recipe" is still worded as if E018 hadn't happened.

## Methodology findings

1. **"does not grow with $n$" (l. 570; caption l. 554 "flat in $n$"; abstract
   ll. 62–63 "stable over the accessible range of $n$"; contribution 4 l. 140
   "scale-stable in $n$") is contradicted by Table II itself.** Regret rises
   monotonically in $n$ in all 7 families at $p{=}1$ (3-reg 0.037→0.055,
   ER(0.25) 0.031→0.045, WS(0.5) 0.031→0.048) and 7/7 at $p{=}3$; with the
   stated SE (0.001–0.006) and the paper's own 2×SE rule (ll. 478–479), the
   $n{=}12\to18$ deltas (up to 0.018) are significant. A 7/7 monotone sign
   pattern is itself $p\approx0.008$. The density law even *predicts* mild
   growth ($m$ grows with $n$ off-regular families). Fix: report the trend
   honestly ("grows mildly, ≈0.002 AR per qubit at $p=1$; no sign of
   compounding") everywhere the claim appears, or add a trend test.
2. **Three inconsistent SE statements.** Setup l. 478: "0.002–0.006"; Table II
   caption l. 554: "~0.001–0.006"; the generated table file's comment: "~0.001".
   Pick one convention and state the max per-cell SE.
3. **Abstract l. 68 and contribution 5 l. 147: sampled $N$ "matches exact
   parameter setting within 0.01"** — the text (ll. 692–694) and E010 say 0.011
   ($p{=}1$) and 0.016 ($p{=}3$). 0.016 is not "within 0.01." Fix: "within
   0.02" or quote both numbers.
4. **"all 28 (family, $n$, $p$) cells" (l. 713).** E018 has 56 cells
   (7×4×2); 28 is the subset where the sampled-$N$ variant exists
   ($n{=}16,18$). Transfer in fact wins all 56 against the exact compression —
   say so: "every cell (56 vs the exact compression; all 28 where sampled $N$
   was measured)". A referee who counts will flag the 28.
5. **The $p{=}3$ "regret roughly doubles" (ll. 568–569) conflates depth with
   ceiling construction.** The $p{=}3$ ceiling is best-of-4096 *linear ramps*
   on an 8-points-per-axis 4-D grid (ll. 470–471) vs 40 per axis in 2-D at
   $p{=}1$; coarser quantization alone inflates both proxy and transfer regret.
   Comparisons are internally consistent (shared grids — good), but add one
   sentence so "doubles" is not read as depth-intrinsic.
6. **Argmax-displacement as "the only predictor" (ll. 651–652, 611–613) omits
   E012's own two caveats**: the correlation is partly mechanical (on a smooth
   landscape the value gap grows with displacement) and the predictor is post
   hoc (it requires the true argmax — a diagnosis, not a certificate). One
   clause fixes it and actually strengthens the thesis ("no *a priori* scalar
   certifies transfer").
7. **"At classically simulable sizes the proxy buys no wall-clock advantage"
   (ll. 755–757; intro ll. 104–106) overstates its own evidence.** E013
   measured $n\le20$ and *extrapolates the pipeline's advantage opening at
   $n\approx22$–$24$* — inside the simulable range. Fix: "at the sizes we
   measured ($n\le20$)" and restore the crossover estimate, which currently
   appears nowhere in the paper.
8. **Timing generality (Table III, ll. 726–729).** One instance per $n$,
   medians of 3; cost scales with $m$. Say "one instance per $n$" in the
   caption. Minor: abstract l. 51 says leakage is "computable in $O(2^n)$";
   the statevector layer it requires is $O(n2^n)$ (text l. 347 is fine).

## Coherence/overselling findings

1. **Contribution 5, l. 148–149: "the recommended recipe wherever
   statevector-cost passes are affordable" is contradicted by the paper's own
   E018 paragraph (ll. 707–713)**: wherever statevector passes are affordable,
   an $n{=}10$ same-family source is trivially evaluable and plain transfer
   "beats every proxy variant … usually by an order of magnitude." Fix: "the
   recommended distribution *for the proxy*, when no same-family source
   exists."
2. **Same tension inside the abstract** (ll. 68–71 vs 71–74): the sampled-$N$
   sentence claims a regime that the very next sentence hands to transfer. The
   reconciling qualifier ("family-agnostic" ⇒ for instances with no usable
   family/source) is implicit; make it explicit, or put the transfer result
   first and scope the recipe after it.
3. **Bolded "Recommended recipe" (ll. 698–699) stands unqualified** one
   paragraph before its defeat. Add "among proxy variants" inside the bold.
4. **Abstract ll. 63–64: "at depth it decouples from state fidelity
   *entirely*"** — measured at $p{=}3$ only ($\rho=-0.02$, ll. 596–599), while
   the proxy's advertised regime is $p{=}20$. Write "already at $p=3$".
5. **The hardware-regime claim (ll. 109–110, 152–153, 719–722, 759–760) is
   asserted, never defended against E018's own logic.** A hardware target of a
   known family still has a cheap classical source landscape (simulate
   $n{=}10$, transfer — zero shots spent on parameter setting). The honest
   statement: at extrapolation ranges far beyond the tested $10\to18$, *both*
   transfer and the analytical proxy are unvalidated; the proxy's residual
   regime is instances with no family structure. Say that, or a referee will.
6. **Nowhere does the paper flag that its parameter-quality evidence stops at
   $p=3$** while §I (l. 90) advertises the heuristic's $p{=}20$ use case and
   leakage is tracked to $p{=}30$. One sentence in §V-A or Limits.
7. Draft artifacts that must not reach submission: "the literature pass
   predates final submission" (l. 779) — an internal TODO in the Limits
   paragraph — and the header "not for distribution" / "Author list to be
   finalized" (ll. 36–38).

## Cut/structure suggestions

- **Abstract: ~330 words → ≤200.** It is currently a compressed results
  section (a dozen numbers, em-dash triple-clauses, e.g. ll. 64–67). Keep:
  exactness, error calculus, third-order law, density, argmax transfer, the
  transfer-beats-proxy honesty. Numbers can go.
- **Split §V-D (ll. 616–722), which does four jobs** (analytical robustness +
  filters; fitted-shape dissection E012+E017; sampled-$N$ recipe; E018
  transfer). Move recipe + E018 + §V-E timing into one closing "Practical
  accounting" subsection ending with the regime map — the recipe and its
  defeat then sit honestly side by side instead of whiplashing.
- **Compress Lemmas 4–5 + Remark (ll. 435–460)** to a 6-line remark citing the
  repo (the $\tau^2/m$, Pearson-0.994 detail is repo material). Buys ~0.3 pp
  for the fixes above.
- **Merge the trajectory-PCA paragraph (ll. 539–547) with the E016 Discussion
  paragraph (ll. 783–793)** — two places currently tell one "mis-aim /
  rotating subspace" story; keep the Discussion version short.
- **Add a 4-sentence Conclusion.** The paper ends on an open problem; QCE
  readers expect a restated thesis + regime map. The Discussion otherwise does
  real work (E016 is a good measured obstacle) — keep it.
- Dense-sentence pass (old T3.1) is still pending per STATUS; worst remaining
  offenders: ll. 64–67, 681–687, 699–705.

## Top 5 unanswered referee questions (by risk)

1. **Weighted/continuous costs.** For weighted MaxCut (the case of
   \cite{shaydulin2023,sureshbabu2024}, which the paper cites) cost classes
   are generically singletons: $\Hhom$ is the full space, the compression is
   trivial, and the recursion is $2^n$-dimensional. Does anything survive
   binning, and what does binning do to Thm 1's exactness? Not mentioned once.
2. **What is the deliverable in the advertised regime?** Beyond simulation,
   every new instrument here (leakage, sampled $N$, empirical $P$, regret
   measurement) is unavailable; what remains is sud2024's analytical $N$ plus
   this paper's warning — while transfer from a small simulated source also
   remains available. Concretely, what should a practitioner at $n{=}100$ do
   differently after reading this paper?
3. **Depth gap.** Parameter-quality claims stop at $p{=}3$ ramps; the proxy's
   selling point is $p{=}20$. Does argmax transfer persist, and does regret
   keep "doubling," at $p{=}10$–$30$? (Leakage to $p{=}30$ exists; regret does
   not.)
4. **Noise.** The hardware regime is claimed but all dynamics are exact and
   unitary; the objective on hardware is shot-noisy and biased. Does argmax
   robustness survive realistic sampling noise / depolarization, even in a
   toy model?
5. **Beyond MaxCut.** Lemma 1 (the cancellation mechanism) is MaxCut-specific;
   sud2024 covers Max-$k$-XOR and SAT. Does first-order leakage cancellation
   generalize to other $k$-local costs, and do the density law and argmax
   findings hold on any non-MaxCut problem? A single small experiment would
   defuse this.
