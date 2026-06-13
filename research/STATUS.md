# Research status

*Last updated: 2026-06-13 (overnight autonomous session). This is the one page
Spencer needs to read. Everything here links to a reproducible experiment or a
committed document.*

## Overnight digest (2026-06-13)

Claude worked the paper autonomously while you slept. Everything below is on the
`ClaudeResearch` branch of each repo, **committed locally, not pushed, paper not
merged to main**. Highlights:

1. **`research/paper_explainer.md`** — plain-language companion (you said the paper
   was hard to follow): the proxy-as-compression reframe, the two error axes, the
   three theorems as a story. Read this first.
2. **`research/paper_critique.md`** — adversarial review + a full number audit
   (every quoted statistic checks out; two errors found and fixed).
3. **§5.3 now has its missing figure + regret table** — the central claim was
   prose-only. Figure shows the leakage→regret ranking honestly: ρ=0.96/0.86 at
   p=1, ≈0 at p=3. Paper compiles clean (14 pp). [exp 004]
4. **Three unambiguous draft fixes applied** (−0.007 instance mis-attributed;
   E001–E012→E013; a depth-overlap phrasing).
5. **New result E015** — opened the V₂ black box behind the density law: two exact
   lemmas (verified 1e-10) + a triangle-conditioning mechanism. Paper-ready LaTeX in
   `research/v2_density_law.md`.

**"Decisions needed from Spencer" (bottom) is no longer empty** — five items shape
the paper; none block further work.

**Goal:** QCE 2027 contributed paper, *"When and why does the homogeneous proxy work?
QAOA parameter setting as subspace compression."* Complete draft by **Aug 15, 2026**;
arXiv late August; submit April 2027. Full plan: see `research/decisions.md` and the
program plan (Claude's plan file, to be mirrored into `research/program.md`).

## Where we are

**Phases 1–3 complete; Phase 4 (assembly) underway.** The paper repo
(`papers/OverleafPaper`, branch `ClaudeResearch`) now holds
`theory_compression.tex` (standalone theory notes with full proofs) and
`qce2027_paper.tex` (full IEEE-format draft: theory + experimental anatomy
from experiments 001–013 + figures; 4 pp with room to grow). E2.4 is now done
(exp 012, below) and folded into the paper's model-error section. The
systematic literature pass is done (`research/literature_pass.md`): **all four
claimed-new items survive** — no overlapping prior work found among Sud et
al.'s ~22 citers or the projection/lumping/MOR/surrogate literature; five
citations added to the paper (closest neighbor: Krüger–Mauerer landscape
approximation, Quantum 2025 — complementary, cited). One re-check recommended
at submission time. Remaining: figure/table expansion, author list and
acknowledgments. Overleaf shows only `main`, so
merge the paper repo's ClaudeResearch branch when ready to edit there.

## What we know (established results only)

- **The homogeneous proxy is exactly an orthogonal subspace compression of QAOA**
  (not an approximation): running the proxy with the same-instance empirical
  N(c';d,c) equals evolve-one-layer-then-project-onto-cost-classes, to machine
  precision. → [experiments/001_proxy-is-compression](experiments/001_proxy-is-compression/README.md)
  Consequence: total proxy error splits cleanly into *compression error*
  (leakage out of the cost-class subspace, model-independent) + *model error*
  (using analytical/fitted N instead of the instance's own).
- **Parameter setting is well-posed, and the compression itself is nearly lossless
  at p=1 on every graph family** (ER, BA, WS, 3-regular; n=12–14, 420 instances):
  value-added over the balanced-partition baseline is positive for every single
  instance, and the exact-compression proxy lands within ~0.03 AR of the true grid
  ceiling on all families. Regret grows to ~0.06–0.09 at p=3 ramps.
  → [experiments/002_baselines-and-headroom](experiments/002_baselines-and-headroom/README.md)
  Two consequences: (i) the old log's "proxies fail on non-ER graphs" must be a
  *model-error* story (bad analytical N), not a compression-error story, at least
  at p=1 — sharpens what E1.3 must separate; (ii) depth, not graph family, is
  where compression error first bites.
- **The Theorem-2 bound is tight enough to use, and leakage predicts fidelity.**
  At p=20 ramps (840 runs, n=12–16): ‖ψ−φ‖ ≤ Σλ holds with median slack ≈ 4×
  (never >10×), Σλ is a near-functional predictor of the actual distance, and it
  ranks families by fidelity at Spearman ρ ≈ 0.96–1.0 (except the decayed
  n=16/large-angle cell). → [experiments/003_leakage-vs-overlap](experiments/003_leakage-vs-overlap/README.md)
- **Surprise: density, not "ER-ness," drives compression error** — family-mean Σλ
  tracks edge count at Pearson 0.97; ER(0.5) is the *worst*-compressing family
  tested, 3-regular the best (at fixed unscaled angles; see exp 003 caveats).
- **GATE (exp 004): leakage ranks families by parameter-setting regret at p=1**
  (Spearman ρ = 0.86–0.96, criterion ≳0.8 met); at p=3 along proxy-chosen ramps
  both quantities barely vary across families (ρ ≈ 0 — no signal, not a
  contradiction). → [experiments/004_gate-leakage-vs-regret](experiments/004_gate-leakage-vs-regret/README.md)
- **H-ER-specific (old log) is contradicted**: analytical PaperProxy with effective
  edge probability gets regret ≈ 0.01–0.02 at p=1 on BA/WS/3-regular/sparse-ER —
  *better* than exact compression (the smooth model regularizes parameter choice).
  Its one failure is dense ER(0.5), where its argmax lands on an *unphysical
  norm-inflated artifact* (predicted ⟨C⟩ = 93 on a 38-edge graph). Since exact
  compression is contractive (Thm 1), norm inflation certifies model error for
  free; a norm sanity filter repairs the diagnosed instance
  (`004/diagnose_paper_artifact.jl`); exp 005 (running) quantifies strict-vs-loose
  thresholds across all instances.
- **Timing (exp 013): the exact N is the pipeline's wall, and brute force is
  free at simulable sizes.** On an A100, exact N costs 21.8 s at n=20 (clean
  O(4ⁿ) growth) vs 0.48 s for the brute-force p=1 ceiling it competes with;
  sampled N (S=10) is 50× cheaper (0.44 s). GPU statevector passes are
  kernel-launch-bound (~0.2–0.3 ms, flat in n ≤ 20), so exhaustive grid search
  matches the proxy pipeline end to end even at p=20; the pipeline's cost
  advantage opens around n ≈ 22–24 (extrapolated). Running n=20 exact N
  required grid-striding the GPU homodist kernel (commit `ce337950`). The
  paper prints the measured table in §5.3.
  → [experiments/013_timing-benchmark](experiments/013_timing-benchmark/README.md)
- **V₂ (the cubic-leakage variance functional) has an exact anatomy (exp 015).**
  Two machine-verified lemmas (1e-10, 280 instances): (L1) V₂ is the within-cost-class
  variance of T(y)=Σ_i(Σ_{j∼i}s_j)², and (L2) its *unconditional* variance is exactly
  2·Σ_{j≠k}A_{jk}² (squared codegrees), → 8p²m² for ER. The cubic law λ₁=(βγ²/8)√V₂
  holds to <0.3% everywhere. The density law √V₂∝m is an approximate cancellation:
  codegree variance grows super-linearly with density, conditioning on cost removes a
  triangle-driven fraction (ρ≈0.37 dense ER → 0.80 sparse 3-regular), and the two
  nearly cancel. Advances the open V₂ problem to one quantity (the conditioning
  correction). → [experiments/015_v2-density-law](experiments/015_v2-density-law/README.md);
  derivation + paper-ready LaTeX in [research/v2_density_law.md](v2_density_law.md)

## Working hypotheses (NOT established — from the deleted research log or intuition)

- H1: The proxy's usefulness is governed by leakage out of the cost-class subspace,
  which concentrates only for ER-like graphs. (The paper's central claim.)
- ~~H2: Fitted proxy shapes (triangle/Gaussian) with lower entrywise MSE against the
  empirical N(c';d,c) gave *worse* parameter setting.~~ **Confirmed** (exp 012:
  Triangle fitting raises regret on 136/140 instances; Normal fitting improves the
  objective 10× and never changes the chosen parameters). The proposed
  dynamics-weighted *amplitude* norm does NOT explain it — regret is an
  argmax-transfer quantity (see exp 012).
- H3: Sampling 5–10 bitstrings per cost class suffices to estimate the quantities that
  matter. (Old log; needs re-verification in the leakage metric, not entrywise N.)
- ~~H4: On ER(0.5), random balanced partitions reach ~75–85% approximation ratio.~~
  **Confirmed** (exp 002: 0.75–0.76 mean on ER(0.5) and BA(k=4) at n=12–14); folded
  into metric design as "value-added" = AR_proxy − AR_baseline.
- H2/H-ER-specific is now under pressure: exp 002 shows the *compression* is
  near-lossless on non-ER families at p=1, so any non-ER failure must come from the
  analytical/fitted N (model error). E1.3 tests this directly.

## Open questions / next experiments

- ~~E0.1: is the proxy numerically identical to the compressed evolution?~~ **Done, yes** (exp 001).
- ~~E1.1: is there headroom above baselines?~~ **Done, yes, on every family** (exp 002).
- ~~E1.2: does Σλ track the overlap deficit?~~ **Done, yes — bound slack ~4×,
  near-functional predictor** (exp 003).
- ~~E1.3 gate~~ **Run: PASS at p=1 (ρ = 0.86–0.96); p=3 inconclusive** (exp 004).
- ~~E005: norm filter as universal recipe~~ **Done: NEGATIVE — rescues only dense
  ER(0.5) (regret 0.159 → 0.077) and is catastrophic elsewhere; the analytical
  model's calibration breaks long before its argmax moves.** The argmax-robustness
  of the analytical proxy is itself a finding for §6.
  → [experiments/005_norm-filtered-paper-proxy](experiments/005_norm-filtered-paper-proxy/README.md)
- ~~E006: physicality cap~~ **Done: NEGATIVE — the model inflates predictions even
  at its (correct) peak on sparse families, so any value-based veto rejects the
  answer; on dense ER(0.5) the spurious region extends below the physical cap.**
  Filter arc (004→005→006) closed: the analytical proxy's values are pure noise in
  absolute terms; its argmax location is the only usable signal — excellent off
  dense graphs, corrupted on dense ER(0.5) by a spurious large-β peak.
  → [experiments/006_physicality-filter](experiments/006_physicality-filter/README.md)
- ~~E007 = E2.1: leakage anatomy~~ **Done — three results:** (i) small-angle law
  λ ≈ const·β·γ²·m, with the O(βγ) term killed by an exact MaxCut identity
  (Σ_i c(x⊕e_i) = (n−4)c(x) + 2m — now a committed test) → quantitative
  explanation of the proxy's small-γ accuracy, lemma for §3/§4; (ii) density law
  sharpened: λ/(βγ²m) constant to ±40% across families; (iii) compression leakage
  does NOT flag the analytical proxy's spurious peak (38th percentile) —
  compression and model error are empirically independent axes, a central claim.
  → [experiments/007_leakage-anatomy](experiments/007_leakage-anatomy/README.md)
- ~~E2.2 sampled-leakage predictor~~ **Done: S=5 estimates η_F to median 3.2%
  error (uniform across families/angles) — far below the between-family
  differences; Phase-3 scale-up enabler confirmed. Bonus: Theorem 3's identity
  verified to 4.5×10⁻¹⁵ over 210 full-enumeration rows.**
  → [experiments/008_sampled-leakage-predictor](experiments/008_sampled-leakage-predictor/README.md)
- ~~E2.3 trajectory-PCA~~ **Done: "wrong subspace," decisively — p=20 ramp
  trajectories have effective dimension 2–4 (99% energy), so a 2–3-dim PCA
  subspace matches the entire (m+1)-dim cost-class subspace; proxy degradation
  is mis-aim of the fixed cost-class frame, not state complexity.**
  → [experiments/009_trajectory-pca](experiments/009_trajectory-pca/README.md)
- **E3.1 (scale-up, n=16–18): the structure survives scale unchanged** — regret
  flat in n (p=1: 0.03–0.05; p=3: 0.08–0.11), value-added positive everywhere at
  family level, and **the sampled-N proxy (S=10) matches exact-N parameter
  choices within ~0.01 AR with slightly LOWER regret** — the practical recipe
  the filter arc failed to find (sampled N + empirical P). Family regret
  differences compress at scale, so ranking claims must be scoped.
  → [experiments/010_scaleup-ranking](experiments/010_scaleup-ranking/README.md)
- **E3.2 (depth, p=30, n=16–20): leakage accumulates linearly in depth,
  quantitatively as Theorem 3 predicts** (ratio 1.505 vs 1.5 at small ramps),
  sublinearly in m along deep trajectories, with p=30 small-ramp overlaps of
  0.71–0.81 even at n=20.
  → [experiments/011_depth-scaling](experiments/011_depth-scaling/README.md)
- ~~E2.4 fitted-shape paradox~~ **Done: H2 confirmed strongly, and the planned
  explanation is refuted in its simple form.** Entrywise-MSE fitting of shapes is
  harmful (Triangle: mean p=1 regret 0.054 → 0.211, worse on 136/140) or inert
  (Normal: 10× better MSE, argmax identical to the unfitted default on 140/140).
  No scalar mismatch norm tested predicts regret (entrywise MSE ρ≈−0.1,
  one-layer amplitude error ρ≈−0.2, landscape Pearson ρ≈+0.1); argmax
  displacement does (ρ≈0.7). Raw PaperProxy has amplitude error ~1e5 (slice sums
  up to 9e8) yet near-best regret — amplitude norms are hostage to scale
  conventions the argmax ignores. Per-slice calibration to 2^n is argmax-neutral
  at p=1 (0/140 changes for shapes). Extends the 004→005→006 arc: for *every*
  model class, argmax location is the only usable signal; §6 recommendation is
  "don't fit shapes by MSE — use defaults or analytical N."
  → [experiments/012_fitted-shape-paradox](experiments/012_fitted-shape-paradox/README.md)
- Next up: T2.0 (LaTeX write-up of Thms 1–3 + the O(βγ²) lemma); systematic
  lit pass (in progress, `research/literature_pass.md`). Then Phase-4 assembly.
- Open framing question for the paper (not blocking): at p=3 the proxy-chosen
  schedules equalize leakage across families — is regret there governed by
  landscape/argmax robustness rather than fidelity? Phase-2 experiments decide
  how §5 narrates depth.

## Decisions needed from Spencer

Surfaced in the overnight critique (`research/paper_critique.md`); none block
ongoing work, but they shape the paper:

1. **Author list + acknowledgments** (paper lines 39, 802) — only you can fill these.
2. **Venue/format**: the draft header says "IEEE TQE" (a journal) but the plan says
   QCE 2027 (a conference: IEEEtran, ~8–12 pp, different review). Which target?
   It gates reformatting and the length budget.
3. **Headline framing**: the leakage→regret ranking is the paper's *weakest*
   evidence (7 points, p=1-only — collapses to ρ≈0 at p=3). Critique T1.1 recommends
   leading with the mechanism + argmax-transfer and demoting the ranking to
   corroboration. Agree? (The new §5.3 figure already scopes it honestly.)
4. **Insert the E015 V₂ lemmas into §4?** They are exact and machine-verified; a
   ready LaTeX block is in `research/v2_density_law.md`. Left §4 untouched pending
   your OK — the theory is yours to own.
5. **`theory_compression.tex` fate**: retire it, or keep as the "full proofs"
   companion/appendix? It now duplicates §3–4 and can drift out of sync.
