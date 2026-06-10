# Research status

*Last updated: 2026-06-10 (evening). This is the one page Spencer needs to read.
Everything here links to a reproducible experiment or a committed document.*

**Goal:** QCE 2027 contributed paper, *"When and why does the homogeneous proxy work?
QAOA parameter setting as subspace compression."* Complete draft by **Aug 15, 2026**;
arXiv late August; submit April 2027. Full plan: see `research/decisions.md` and the
program plan (Claude's plan file, to be mirrored into `research/program.md`).

## Where we are

**Phase 1 gate has run (planned ~Jul 8, ran Jun 10): PASS at p=1, inconclusive
at p=3** — the theory-led framing proceeds, with the depth question moved to
Phase 2's schedule-diverse measurements. One follow-up sweep (norm-filtered
PaperProxy, exp 005) running in background.

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

## Working hypotheses (NOT established — from the deleted research log or intuition)

- H1: The proxy's usefulness is governed by leakage out of the cost-class subspace,
  which concentrates only for ER-like graphs. (The paper's central claim.)
- H2: Fitted proxy shapes (triangle/Gaussian) with lower entrywise MSE against the
  empirical N(c';d,c) gave *worse* parameter setting. (Old log; needs re-verification.)
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
- Next up: T2.0 (LaTeX write-up of Thms 1–3 + the new O(βγ²) lemma + lit pass,
  per `research/theory_notes.md`); E2.2 sampled-variance predictor; E2.3
  trajectory-PCA subspace comparison. Phase 2 mechanism work, then scale-up.
- E2.1 (next up): (γ, β) leakage anatomy per family — separates the density effect
  from angle rescaling (H-density), tests Theorem 3's variance identity, and
  supplies the schedule diversity the p=3 gate lacked.
- Open framing question for the paper (not blocking): at p=3 the proxy-chosen
  schedules equalize leakage across families — is regret there governed by
  landscape/argmax robustness rather than fidelity? Phase-2 experiments decide
  how §5 narrates depth.

## Decisions needed from Spencer

*(none right now)*
