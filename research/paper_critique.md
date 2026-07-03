# Adversarial critique of the QCE 2027 draft

*Written for Spencer, 2026-06-13, reading `qce2027_paper.tex` and
`theory_compression.tex` as an unsympathetic-but-fair QCE/IEEE reviewer would.
Prioritized and actionable. A separate background number-audit cross-checks every
quoted statistic against the experiment records; its results are folded into
§5 below.*

## Top-line assessment

This is a genuinely good paper with a real idea. The reframe (proxy = exact
orthogonal compression of QAOA) is correct, clarifying, and — as far as the
literature pass found — new. The error calculus is elegant and *cheap to
measure*, the MaxCut cancellation lemma is a clean quotable result, and the
honest negative results (filters can't repair model error; MSE-fitting
backfires) are the kind of thing reviewers respect and competitors omit. The
reproducibility scaffolding is better than 95% of QC papers.

**The danger is not the science; it's the framing and the presentation.** Two
things could draw a weak-accept-to-reject from a tough reviewer: (1) the
*headline* empirical claim (leakage ranks families by regret) is the paper's
*weakest* evidence — a 7-point rank correlation that operates at the edge of
statistical resolution and only at p=1 — yet it is sold as central; and (2) the
prose is dense enough that the reviewer has to work to find the contributions
(this is the same wall Spencer hit). Both are fixable without new experiments.

Below: Tier 1 (fix before submission), Tier 2 (strengthen the argument), Tier 3
(polish), Strengths to preserve, and the number audit.

---

## Tier 1 — fix before submission

### T1.1 — The headline claim is the weakest evidence. Reframe what's "central."

The abstract and intro lead with the family-ranking correlation
(leakage ⇔ fidelity ⇔ regret). But by the paper's own §5.3:

- It is a **Spearman over 7 points** (7 family means). ρ = 0.96 (n=12), 0.86
  (n=14). Over 7 ranks, ρ=0.86 is suggestive, not decisive (two-sided
  p≈0.02), and it's a correlation of *family means*, so the effective sample is
  structurally 7 no matter how many instances feed each mean.
- It holds **only at p=1 and small n**; the paper states the signal "fades" at
  deeper p and larger n as family regret differences compress to a few 0.01.
- The family differences at p=1 (spread ~0.009 AR) sit ~9× above the actual
  per-cell SE (~0.001), so the *extreme* families (ER(0.5) lowest, 3-regular
  highest) are cleanly ordered. But the **middle five families are within 1–2 SE
  of each other**, so most of the 7 ranks are not individually resolvable —
  which is exactly why ρ is 0.86–0.96 and not 1, and why it is fragile.
- **The ranking is p=1-only.** Recomputed from exp 004/010: ρ collapses to
  **0.04 (n=12) and −0.04 (n=14) at p=3**, and is only 0.61/0.29 at p=3 n=16/18.
  The signal does not merely "fade" with depth — at p=3 small n it is *gone*.

A reviewer will circle this and write "the central empirical claim is a 7-point
correlation near the resolution limit, valid only in the easiest regime." That's
avoidable, because **it isn't actually the paper's strongest result.** The strong
results are: the *exact* reframe + error calculus (Thms 1–3), the *proved*
third-order MaxCut cancellation, and the *argmax-transfer* reorganization of the
whole parameter-setting question (which makes the negative results coherent).

**Recommendation.** Demote the family-ranking from "central claim" to "one
consistent line of evidence, resolvable at p=1." Let the title question — *when
and why does the proxy work?* — be answered primarily by the **mechanism**
(leakage = measurable compression error; third-order protection; density sets
the constant) and by **argmax-transfer**, with the ranking as corroboration.
This is a re-emphasis, not a retreat: every sentence stays true, but the paper
stops staking its headline on its thinnest plot. (See T1.2 for the figure that
makes the scoped version honest.)

### T1.2 — §5.3 has no figure and no table. The thesis paragraph is bare prose.

§5.3 ("Compression error and parameter quality") is where the leakage→regret
claim lives, and it is **pure prose** — two Spearman numbers in a sentence —
while four *supporting* results each get a figure. The program plan called for a
"headline ranking figure + regret table" here.

The tension: a 7-point scatter looks thin (which is presumably why it's prose),
but leaving the thesis figureless makes it look unsupported. Resolve it by
showing the claim **with its honest scope built in**:

- A **regret table** (family × n × {p=1, p=3}) with instance standard errors —
  turns a paragraph of inline numbers into something scannable and lets the
  reader see the family spread *and* its overlap with the error bars. Low risk,
  pure win. Data: exp 002 (n=12–14) + exp 010 (n=16–18).
- A **ranking figure**: leakage-rank vs regret-rank for the 7 families at n=12
  and n=14, with the Spearman ρ annotated, plus a second panel at p=3 / larger n
  showing the differences *collapsing into the error bars* — i.e. plot the fade
  honestly rather than describing it. That converts the weakness into a
  controlled, scoped statement and pre-empts the reviewer. Data: exp 004 + 010.

(I am building both tonight as task #4.)

### T1.3 — Resolve the venue/format ambiguity; it changes the length budget.

The file header says *"target arXiv then IEEE TQE"* (a journal). The research
plan (`research/program.md`, memory) says **QCE 2027 contributed track** (an IEEE
*conference*, ~8–12 pp, double-column IEEEtran, April 2027). These are different
animals: different template, length, and review culture. The current draft is
`article` 11pt single-column with 1.1in margins — neither target's format.

This is a **decision for Spencer**, and it gates real work (reformatting to
IEEEtran, cutting/expanding to the page limit, possibly splitting theory into an
appendix). Flagging now because it's cheap to decide and expensive to defer.

### T1.4 — "E001–E012" vs the existing, used E013. Consistency.

The abstract, contributions item 5, and "Code and data availability" say
experiments **E001–E012**, but **E013** (timing) exists and *is* the source of
Table 2 (`% E013` comments confirm). Update the ranges to E001–E013 everywhere,
and make sure the code-availability paragraph promises a script for every figure
*and* table actually in the paper (including the new §5.3 ones from T1.2).

---

## Tier 2 — strengthen the argument

### T2.1 — The "density, not ER-ness" headline has a confound the paper under-weights.

λ/(βγ²m) ≈ const is measured **at fixed angles**, and the paper notes (lines
~556–561) that λ ∝ γ²m is invariant under γ↦γ/√m — i.e. the rescaling
practitioners actually use *absorbs the density effect*. So a hostile reading is:
"your density law is an artifact of holding γ fixed; under the schedules people
run, it disappears." The caveat is present but buried mid-paragraph after the
bold claim. **Fix:** state the fixed-angle scope in the same breath as the claim,
and consider a small panel showing leakage under the √m-rescaled schedule so the
reader sees exactly what survives. The honest version ("density drives leakage at
matched angles; matched *schedules* partly cancel it") is still interesting and
much harder to attack.

### T2.2 — "The wrong model beats the exact one" needs a mechanism, not just "regularizes."

§5.4's most striking claim: the analytical N gets *lower* regret (0.01–0.02) than
the exact-empirical compression (0.03–0.06) off-ER, because it "regularizes
instance noise." This is counterintuitive enough that a reviewer will demand
*why* a provably-wrong N chooses better angles than the exact one. The likely
mechanism — the single-instance empirical landscape is ruggeder, so its argmax is
a noisier estimator of the "good angle" than the smooth class-level landscape —
is testable (landscape ruggedness / argmax stability vs smoothing). One sentence
of mechanism + one supporting number would turn a surprising anecdote into a
result. This is exactly the open framing question in STATUS and the seed of my
planned E014.

### T2.3 — Pre-empt the "so the proxy buys nothing below n≈22" reading of Table 2.

The timing table is admirably honest, but a hostile reviewer reads it as: *at
every simulable size, exhaustive grid search matches the proxy pipeline
end-to-end; the cost advantage only opens at n≈22–24.* The paper says the value
at these sizes is "the error calculus and the reuse of one N across depths/grids"
— but that rebuttal is one clause late in §5.4. **Move the framing forward:** in
the intro, set the expectation that the contribution is *understanding* (when/why
it works) plus the high-depth/large-n and hardware regimes, **not** a laptop
speedup. Then Table 2 reads as honesty, not self-sabotage. (Right now the abstract
oversells "enabling classical parameter setting at depths where direct simulation
is intractable," which Table 2 partly undercuts at the tested sizes.)

### T2.4 — Give the reader one parameter-transfer number for calibration.

The paper deliberately uses grid ceilings rather than parameter transfer as the
reference ("so our regrets are conservative"). Defensible. But transfer is *the*
competing cheap method, and a reviewer wants to know how the proxy's regret
compares to just transferring angles from a small instance. The old final_report
has transfer experiments. Even one calibration row ("transfer regret X vs proxy
regret Y on the same cells") would close a predictable referee question. Optional
but high-value.

### T2.5 — The linear-in-depth law: be explicit about what's proved vs fit.

The accumulation ratio 1.505≈1.5 is presented as "predicted." But the prediction
rests on "λ_ℓ depends on (γ_ℓ,β_ℓ)≈f(ℓ/p) alone" — exact only at layer 1 (input
|+⟩); at layer ℓ the input is φ_{ℓ-1}, carrying the whole history. So the
linear-in-p law is a *small-leakage approximation validated empirically*, not a
corollary of the cubic lemma. Say so plainly (one sentence) — otherwise a careful
reviewer catches the gap between "Corollary" and "measured 1.505" and distrusts
the rest.

---

## Tier 3 — polish (also addresses "hard to understand")

### T3.1 — Break the dense sentences. This is a reviewer issue, not just Spencer's.

The abstract's final sentence, and most of §5.4, pack 3 clauses + 2 parentheticals
each. Reviewers skim; buried contributions read as absent. Concretely: split the
abstract's last sentence; give §5.4 a one-line thesis up front ("Across every
model class we tested, only the argmax transfers; no value- or norm-based quantity
certifies it"), then march through the evidence; move inline number-dumps into the
regret table (T1.2) so the prose can *narrate* rather than *list*.

### T3.2 — Decide the fate of `theory_compression.tex`.

It is a standalone doc with its **own** bibliography and **own** numerical-
verification table, overlapping qce2027_paper.tex §3–4. Its header says "intended
to become Section 3," but §3–4 already exist in the main file. Right now it's a
divergence risk (two copies of the theorems that can drift; e.g. confirm the
slack "3–6×" and exponents match between the two). **Decide:** retire it, or keep
it explicitly as the "full proofs" companion/appendix and add a one-line pointer
from the main paper. Don't leave two un-synced sources of the same theorems.

### T3.3 — Tighten the cubic-corollary proof sketch.

The sketch is correct but terse. It would be airtight (and more convincing) if it
noted the load-bearing fact that **|+⟩^⊗n is an exact eigenstate of B(β)** (so
there is zero leakage at γ=0 for *any* β — forcing leakage to carry γ², not just
β). That one observation makes "the first survivor is k=2" obviously complete
rather than asserted, and rules out a competing O(β²γ⁰) term by inspection.

### T3.4 — Make the proved/measured line unmissable in §4.3.

§4.3 ("why homogeneity happens") is the best-written section, but its momentum can
make the *measured* rung (V₂ self-averaging) read as nearly proved. The text does
flag it (~lines 494–497); just sharpen to one explicit sentence: "Rungs 1–4 are
theorems on any graph; rung 5 (self-averaging of V₂) is the one empirical
assumption, and the open problem of §6."

### T3.5 — Glossary is excellent; cash it in earlier.

Table I is genuinely the best on-ramp in the paper. Reference it from the abstract
or first paragraph ("notation in Table I"), and consider pulling its three
"in words" gems (proxy = QAOA + projection; regret; argmax transfer) into the
intro so the reader has the mental model *before* §3.

---

## Strengths to preserve (don't "fix" these)

- The **reframe** and the **two-error decomposition** — the whole contribution.
  Keep them front and center.
- The **cheapness** of leakage (O(2ⁿ), one statevector layer) — emphasize; it's
  what makes the calculus a tool, not a tautology.
- The **negative results** (filters fail; MSE-fitting backfires; no value-norm
  certifies argmax). Rare and credible. Keep reporting them in full.
- The **argmax-transfer** reorganization — it's the idea that makes §5.4 cohere.
- The **reproducibility** scaffolding and per-experiment READMEs.
- The **Lemma** (neighbor-cost sum) — a clean, citable, graph-agnostic identity.

---

## Number audit (paper claims vs experiment records)

A background agent cross-checked **every** quoted statistic against the
experiment READMEs and recomputed headline numbers from the `results*.csv`. The
verdict is reassuring: **the paper is highly accurate.** ~20 claims checked;
Tables 1 and 2 reproduce *exactly*, the Spearman ρ (0.964/0.857), the bound
slack (max 5.80×), the accumulation ratio (1.505), the variance-identity
tolerance (4.48e-15), the cubic exponents (0.978/1.95 vs stated 0.98/1.92), the
fitted-shape numbers, the sampled-N recipe, and the timing table all match the
recorded data. No overreach was found between the experiment READMEs' cautious
"Answer" lines and the paper's statements — if anything the paper is careful.

**Two real fixes (both no-input-needed, slated for task #5):**

1. **Factual error — the "−0.007" instance is misattributed** (lines ~614–616).
   The paper says "one $G(18,0.5)$ instance posts a marginal $-0.007$" inside a
   sentence about the **exact compression**. The data (exp 010) say: it is at
   **$n=16$, not $n=18$**, and its value-added is negative only for the
   **sampled-$N$ proxy** — the *exact-compression* value-added on that instance
   is $+0.013$ (positive). As written the sentence is wrong on both the size and
   the proxy. Fix: either drop it, or move it to the sampled-$N$ paragraph as
   "one ER(0.5) $n{=}16$ instance under sampled $N$."

2. **"E001–E012" → "E001–E013"** (lines 153 and 794). E013 (timing) exists, is
   cited inline (lines 741, 749), and *entirely backs Table 2*. The experiment
   range in the contributions list and the code-availability section must
   include it.

**One phrasing nit (no number error):** "small-ramp overlaps remain 0.71–0.81 at
$p=30$, $n=20$" (line ~584) — 0.71–0.81 is the range across $n=16$–$20$; at
$n=20$ specifically it is 0.71. Rephrase to "0.71–0.81 across $n=16$–$20$ (0.71
at $n=20$)."

**Robustness notes:** the "79–99%" cost-class capture (claim 17) and the γ
exponent 1.92 (vs my 1.95) are both defensible — the former is the README's
family-mean convention (per-instance min dips to 0.74), the latter is sensitive
to the small-angle fit window. Neither needs changing; just be aware if a
reviewer recomputes.
