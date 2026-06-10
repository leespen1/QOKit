# Research Revamp: QCE 2027 Paper + Automated AI Research Workflow

## Context

The G-RIPS 2024 internship produced the triangle proxy for QAOA parameter setting, but
the original speedup pitch is weakened (proxy arrays are precomputed once and amortized),
and the team doubts the homogeneous heuristic itself (seems ER-specific; random balanced
partitions already hit 75–85% AR on ER(0.5)). Spencer wants a paper acceptable to the
**IEEE Quantum Week (QCE) 2027 contributed track** (deadline ~April 2027; QCE 2026 papers
closed April 27, 2026), satisfying the Mitsubishi sponsors, with Spencer contributing
1–2 hr/day of discussion while Claude does most of the execution.

**Decisions made by Spencer (2026-06-10):**
- **Fresh start** on prior experiments. The existing `research_log/` (staged for deletion)
  was produced by an earlier Claude workflow that generated results Spencer didn't
  understand. Its findings are *private hypotheses to re-test*, never citable. The #1
  workflow requirement is **legibility to Spencer**.
- **Venue:** QCE 2027 contributed paper only (no 2026 poster/workshop rush).
- **Framing:** theory-led — *"when and why does the homogeneous proxy work?"*
- **Autonomy:** Claude designs/runs experiments autonomously (automode); Spencer reviews
  digests and makes decisions in daily sessions.
- **Timeline:** complete, submission-ready draft by **mid-August 2026** (~9.5 weeks),
  then arXiv preprint (sponsor deliverable) and QCE 2027 submission in April with ample
  polish buffer. Feasible with trimmed scale-up (see schedule); main risk is the Phase-1
  go/no-go gate (~July 8) forcing a pivot — even then mid-August yields the pivot draft.

## The scientific core (verified against the code)

The homogeneous proxy iteration in `src/QAOA_proxy.jl` is **exactly** the compression
`T = D^{-1/2} (P·U·P) D^{1/2}` of the QAOA layer unitary `U = B(β)Φ(γ)` onto the
(m+1)-dim subspace spanned by normalized cost-class indicator states
`|c⟩ ∝ Σ_{x:c(x)=c}|x⟩`, where `P` is the orthogonal projector and `D = diag(M_c)`
(class sizes) — *when N is the same-instance empirical class average* (what
`get_homogeneous_distribution_from_costs` computes). Analytical/fitted N adds a separate
**model error** on top of this **compression error**. Three provable results anchor the paper:

1. **Theorem 1 (exactness):** proxy = orthogonal compression, with normalization caveats
   (no renormalization in `QAOA_proxy_*`; unattained cost classes; `expectation()` uses
   binomial P which is a third convention — pick one, report discrepancies).
2. **Theorem 2 (telescoping bound):** `‖ψ_p − φ_p‖ ≤ Σ_ℓ λ_ℓ` where
   `λ_ℓ = ‖(I−P)U_ℓ|φ_{ℓ−1}⟩‖` is per-layer **leakage** — cheap to measure (one
   statevector layer + O(2^n) projection; no 2^n×2^n matrices ever).
3. **Theorem 3 (concentration link):** one-layer leakage `λ(c')²` is *exactly* the
   class-size-weighted within-class variance of `g_{c'}(y) = Σ_d f_d(β) n(y;d,c')` —
   i.e. the relevant concentration is of the **f_d(β)-weighted** combination of
   n(x;d,c), not entrywise. Explains the small-γ/β regime and why entrywise-MSE
   shape fitting was the wrong objective. Small-(γ,β) perturbative corollary derivable.

**Central empirical claim to validate:** graph families ranked by leakage ⇔ ranked by
proxy fidelity ⇔ ranked by parameter-setting regret. Sampling-based N estimation and the
fitted-shape negative result become supporting evidence.

**Riskiest assumption (test FIRST, fail-fast):** that state-compression error controls
*parameter-setting quality* — a low-fidelity proxy could still have a good argmax, and
ER(0.5)'s flat landscape (random-partition floor ~75–85% AR) could make rankings noise.
Phase 1 has an explicit go/no-go gate with a pre-committed pivot framing ("leakage
explains fidelity; parameter quality is governed by argmax robustness" — still publishable).

## Research phases (compressed: June 10 → August 15, 2026)

Metrics everywhere: **regret** = AR_gridceiling − AR_proxy and **value-added** =
AR_proxy − AR_randombalanced (never raw AR). Experiments run autonomously every day;
writing overlaps with experiments throughout.

- **Phase 0 (wk 1, Jun 10–17): formalization + sanity.** Scaffold workflow v2; E0.1:
  verify Theorem 1 numerically to 1e-10 (proxy iterate vs statevector-layer-then-project),
  n=8–12 ER, minutes on laptop. Watch the `cis(-γc/2)` and `pi_units` conventions
  (figure4 path). Graph generators + leakage primitives land this week.
- **Phase 1 (wks 2–4, Jun 17–Jul 8): go/no-go.** E1.1 baselines/metric design (n=12–14,
  p∈{1,3}, families: ER(0.5), ER(0.25), Barabási–Albert ×2, Watts–Strogatz ×2, random
  3-regular; 30 instances; true-QAOA 40×40 grids for ceilings). E1.2 leakage-vs-overlap
  bound tightness (extends figure4, p=20 ramps). **E1.3 the central correlation**
  (Spearman of leakage/fidelity/regret rankings; empirical-N proxy isolates compression
  from model error). **Gate ~Jul 8:** proceed theory-led iff family-level |ρ| ≳ 0.8;
  otherwise pivot to the pre-committed fidelity-only framing. T2.0 (Theorems 1–3 in
  LaTeX) starts in parallel this phase — it needs no experimental input.
- **Phase 2 (wks 4–6, Jul 8–22): mechanism.** Verify Theorem 3 identity to machine
  precision. E2.1 leakage anatomy η(γ,β) heatmaps per family (GPU batched, n=14–16).
  E2.2 sampled-variance cheap predictor (S∈{2,5,10,25}; re-verifies the untrusted
  "S=5–10 suffices" hypothesis in the quantity that matters). E2.3 cost-class subspace
  vs optimal trajectory-PCA rank-(m+1) subspace (thin SVD of stacked intermediates —
  discriminates "wrong subspace" vs "not low-rank" on BA/WS). E2.4 (fitted-shape paradox
  via weighted-norm vs entrywise model error) — trim first within Phase 2 if behind.
- **Phase 3 (wks 6–8, Jul 22–Aug 5): scale-up, trimmed.** E3.1 headline ranking figure
  at n=16–18 (sampled N; exact GPU homodist spot-checks; n=20 only if HPC queue allows),
  20–30 instances per family, via Slurm batches modeled on
  `scripts/paper_figures/run_all.sb`. E3.2 depth scaling of Σλ_ℓ for ramps to p=30.
  E3.3 (Chung–Lu heterogeneity knob) is **cut** from the mid-August scope; revisit
  during the Aug–Apr polish window if desired.
- **Phase 4 (wks 8–9.5, Aug 5–15): paper assembly.** Reproducible figure scripts (one
  per figure, DrWatson `produce_or_load`, fixed seeds); robustness passes; complete
  submission-ready draft by **Aug 15**. Then: arXiv preprint late August (sponsor
  deliverable); QCE 2027 submission ~April 2027 after the polish window.

**Feasibility conditions:** daily autonomous runs actually happen; Spencer's 1–2 hr/day
is consistent; MSU HPC access works for the ~weeks-6–8 GPU batches; the Phase-1 gate
passes (a pivot still completes by mid-August, with narrower claims).

**Paper skeleton:** Intro → Background → §3 Proxy-as-compression (Thms 1–3) →
§4 Anatomy of leakage (η(γ,β) heatmaps; λ_ℓ vs overlap; subspace-vs-PCA) →
§5 Leakage predicts heuristic quality (headline ranking fig + regret table) →
§6 Practical consequences (sampled predictor; fitted-shape paradox) → Discussion.

## New code (Julia, all ≲100 lines each, in `src/` with tests)

1. `project_onto_cost_classes(state, costs)` → class coefficients + residual norm (leakage primitive)
2. `layer_leakage(costs, n, homodist, γs, βs)` → per-layer λ_ℓ, norms, overlaps
3. `per_class_leakage(costs, n, γ, β)` (+ GPU batched variant via `gpu_qaoa_statevector_batched`)
4. `proxy_transfer_matrix(N, γ, β)` → explicit (m+1)² T for spectral/model-error analysis
5. `sampled_homogeneous_distribution(costs; S)` + `sampled_class_variance` (stratified per-class)
6. Graph generators: BA / WS / random-regular / Chung–Lu (port from
   `data/data_generation/collect_homodists.py` or use Graphs.jl; update `[compat]`)
7. `trajectory_svd(states)` → captured-energy comparison
8. Experiment drivers under `scripts/experiments/` with DrWatson caching + Slurm wrappers

**Reuse (verified to exist):** `qaoa_statevector(...; return_intermediates)`,
`apply_phase_gate!`/`apply_x_mixer!`, `QAOA_proxy_basic/single/multi`, `expectation`,
`get_homogeneous_distribution_from_costs[_direct]` + GPU variants,
`get_real_distribution_from_costs`, `linear_ramp[_matrix]`,
`erdos_renyi_edges`/`maxcut_optimal` (`scripts/paper_figures/common.jl`),
figure4's `proxy_statevector_from_compressed`, `cpu_multi_proxy_mse`/`gpu_multi_proxy_mse`,
`python/grips/sendai_opt.py:fit_proxy_to_real` (E2.4 only). Python stays a cross-check oracle.

## Workflow v2 (designed around "Spencer must understand everything")

Replace the deleted `research_log/` with a new top-level `research/` directory:

- **`research/STATUS.md`** — one page, always current: phase, what we know (one sentence
  per established result, each linking to its experiment), open questions, and a
  **"Decisions needed from Spencer"** list. This is the only file Spencer *must* read.
- **`research/decisions.md`** — append-only log of Spencer's decisions with dates.
- **`research/experiments/NNN_short-name/README.md`** — one per experiment, in fixed
  plain-language format: **Question** (1 sentence) → **Answer** (1 sentence, bold,
  written only after the run) → Method (a paragraph a non-expert coauthor can follow) →
  Figures inline → Caveats → exact repro command (`julia --project run.jl`, fixed seed).
  Negative results get equal treatment. No result enters STATUS.md without a committed,
  reproducible script.
- **Daily cadence (automode):** each session opens with a ≤15-line **digest** (what ran,
  one-sentence answers, what's blocked on Spencer, proposed next runs). Spencer decides;
  Claude launches long runs in the background during/at end of session (Bash
  `run_in_background`, Slurm for Phase 3) and writes up results next session. On the
  compressed timeline, every day should end with at least one run queued.
- **Standards (fail-fast):** fixed seeds; assertions in scripts (no silent fallbacks);
  every claim traceable to a script; untrusted prior-log findings explicitly labeled
  "hypothesis" until re-verified.
- Save workflow conventions + program summary to Claude's persistent memory so every
  future session starts oriented.

## First implementation session (concrete)

1. Commit the staged `research_log/` deletion (Spencer's decision, already staged) and
   scaffold `research/` (STATUS.md, decisions.md recording today's four decisions).
2. Write `project_onto_cost_classes` + test, then **E0.1** (Theorem-1 exactness check)
   via the Julia MCP session with Revise; log it as `experiments/001_*`.
3. Port BA/WS/random-regular graph generators + tests (needed by all of Phase 1).
4. Draft E1.1 driver script; kick off the first background run.
5. Write memory files (workflow conventions, program map).

## Verification

- E0.1 itself is the verification of the math core (mismatch >1e-10 ⇒ stop, fix conventions).
- New `src/` functions get tests in `test/` (e.g. leakage of a state already in the
  subspace = 0; sampled N → exact N as S→all; BA/WS generators match networkx edge
  counts via PythonCall). Final check with `Pkg.test()` before any PR.
- Every experiment README ends with its repro command; Phase-1 gate results reviewed
  with Spencer before Phase 2 begins.
