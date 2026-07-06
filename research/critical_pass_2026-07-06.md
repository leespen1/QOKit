# Referee report — qce2027_paper.tex (adversarial full read, 2026-07-06)

## 1. Post-surgery coherence (highest priority)
The retraction was executed cleanly. Grepping the whole file for the four
orphan classes:
1. **No surviving "proxy wins/beats/superior" claim.** Every `win`/`beat`/
   `dominat` token is either pro-transfer (L154, L791 "transfer wins every
   cell", L819, L831), an internal proxy-variant comparison (L649 analytical
   `N` "beats the exact compression", L671), a hypothetical (L941–951
   instance-adapted frames "could beat"), or the explicitly-corrected artifact
   (L828 "the proxy appears to win by up to $7\times$" — immediately reversed
   at L833 "matches (not beats)"). None contradicts the new thesis.
2. **No live deleted-niche reference.** `crossover` appears once (L839) *only*
   inside the retraction disclaimer ("mistook the mean-rescaling handicap for a
   genuine crossover; we flag the correction"). No CV², no alpha-sweep, no
   pre-registered boundary as a live claim.
3. **`fig:alpha` fully removed** — zero `\ref`/`\label` hits; all 6 figures +
   3 tables are still referenced, no orphaned floats.
4. **Retracted experiments cleanly dropped from citations.** Cited set is
   E001–E022, E027, E029–E032; E023/E024/E025/E026/E028 (the niche run) are
   gone. **No vestigial experiment supports a dead claim.**

**Orphan count of the four target classes: 0.** Only defect: the reproducibility
footnote L169 says "cited in the text as **E001--E031**", but **E032 is cited**
(L839) and E023–E026/E028 are gaps, so the contiguous-range notation is wrong.
Fix: "E001–E032 (with gaps)" or drop the range.

## 2. Under-claiming / contribution clarity
- Contribution *is* foregrounded: the Contributions list (L116–164) leads with
  Exactness → calculus → small-angle structure → V₂; theory is up front, not
  buried. The "error calculus as instrument" spine keeps this from reading as
  "a useless method." **Not incoherent, not fatally under-claiming.**
- But it **over-hedges to the point of self-harm.** The abstract (L44–74) is
  ~30 lines and its middle is a wall of negatives; the intro pre-emptively
  concedes "the contribution of this paper is **not** a laptop-scale speedup"
  (L108) before any result. Abstract last sentence (L71–73, "value is not speed
  but the error calculus itself…") is a real contribution statement but framed
  by negation — reword to a positive lead ("We give an exact compression
  characterization of the proxy and a cheap, measurable leakage certificate…").
- **Title still fits**: the subtitle "…as Subspace Compression" is the true
  contribution and the paper delivers on "When does it work?" as a
  *characterization*. Mild risk that "Work?" promises a practical regime the
  paper concludes doesn't exist; acceptable given the honest framing.

## 3. Accuracy of the Sud reconciliation (checked vs main.tex L451/L465/abstract)
Substantially faithful. Sud reports per-instance AR差 $-.0037{\pm}.0062$,
$.0164{\pm}.0148$, $.0097{\pm}.0183$ (p=1,2,3), "no statistically significant
differences." Draft L189–194 "statistically indistinguishable… mixed in sign"
✓; monotonic-to-p=20 (L192, Sud L465/abstract) ✓; "no transfer table then
existed" (L194, Sud L465 "we were not able to find previous works…") ✓.
Two nits:
- L190 "**per-instance differences $\le0.02$**" — the ≤0.02 values are the
  *means*; per-instance std reaches 0.018 so individual instances exceed 0.02.
  Say "mean per-instance difference ≤0.02."
- L789 "**in 2022**" conflicts with the bibliography entry sud2024 = *PRR* 6,
  023171 (**2024**) (L1006). It's the arXiv-preprint year; make it
  "when [sud2024] appeared" or "in the 2022 preprint" to avoid a referee flag.
- Minor: L790's "transfer fails at p=20 = absence of a table" slightly
  understates that Sud *did* run a non-proxy outer-loop baseline that "fails to
  converge"; consider "no efficient p=20 ramp-transfer schedule then existed."

## 4. Structural balance & length
- **The paper is over length.** ~54k body chars (L44–985) ≈ 8–9k words + 6
  figures + 3 tables in two-column IEEEtran → comfortably **>10pp**. This is the
  single biggest problem and forces cuts regardless of merit.
- Post-retraction the imbalance is **theory-heavy, not practical-thin**: the
  practical section (L734–852) is still rich (sampled-N, transfer calibration,
  high-depth head-to-head, weighted, refinement). The disproportion is the
  **V₂ quartic anatomy** — Lemmas 8–10 + Prop. condcorr (L449–500) plus the
  open-residue paragraph (L953–957) — deep machinery for a density law already
  shown empirically. **Single highest-value change: move Lemma
  quarticsplit + Prop. condcorr to an appendix (or cut to one sentence),
  reclaiming ~1pp for the abstract/intro trim in §2.** Second: the weighted
  paragraph spends ~13 lines (L822–839) litigating the retracted result; keep
  the correction, compress the post-mortem.

## 5. Verdict
- **QCE (contributed):** *Weak accept* — rigorous, reproducible, genuinely
  novel characterization + honest negative, but must be cut to length and
  reframed positively before camera-ready.
- **PRA:** *Major revision* — the core theorem is an elementary projection
  identity and the headline practical result is now "never beats transfer" on
  MaxCut only; the leakage calculus is real but significance/scope must be
  argued much harder for PRA's bar.
