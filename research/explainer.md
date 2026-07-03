# The research, explained simply

*A plain-language walkthrough of everything we've established so far, kept
current as results land. One section per idea, each linked to the experiment
that proves it. For the one-page program state, see [STATUS.md](STATUS.md).
Last updated: 2026-07-03 (through experiment 013).*

## The one-paragraph version

The homogeneous proxy (from the Sud et al. paper) approximates a QAOA state by
tracking one amplitude per *cost value* instead of one per bitstring. Everyone,
including the original authors, treated this as an uncontrolled approximation
that you justify with plots. Our central discovery is that it is not an
approximation at all — it is an **exact orthogonal projection** (a
"compression") of the true QAOA evolution onto a small subspace, *when you
build it from the right N*. That reframing gives us error bars: the amplitude
that "leaks" out of the subspace each layer is a number we can measure cheaply,
it provably bounds the proxy's total error, and it empirically predicts where
parameter setting works and where it fails. The paper tells that story:
when and why does the proxy work? Answer: exactly when leakage is small, and
we can say precisely what makes leakage small.

## The three theorems, in plain language

**Theorem 1 — the proxy is a projection, not a guess.** Take the true 2^n-dim
QAOA state, evolve it one exact layer, then replace every amplitude by the
average over its cost class ("all bitstrings with cost 17 now share one
amplitude"). That averaging step is an orthogonal projection — the same
operation as dropping the off-axis part of a vector. One proxy step is
*exactly* "evolve, then project," provided N(c';d,c) is the empirical class
average from that same graph. We verified this numerically to machine
precision (~1e-15). → [exp 001](experiments/001_proxy-is-compression/README.md)

Why you should care: every error the proxy makes now splits into two separate,
measurable pieces — **compression error** (the projection throws information
away; a property of the graph, not the model) and **model error** (you used a
formula or a fit for N instead of the graph's own empirical N).

**Theorem 2 — leakage bounds the total error.** Each layer, some amplitude
"leaks" out of the cost-class subspace; call its size λ_ℓ. Adding up the
per-layer leakages bounds how far the proxy state can be from the true state:
‖ψ − φ‖ ≤ Σλ. Measuring λ costs one statevector layer — no giant matrices.
In practice the bound is only ~4× from tight, and Σλ is nearly a *function*
of the true error, so leakage is a working error bar, not just a bound.
→ [exp 003](experiments/003_leakage-vs-overlap/README.md)

**Theorem 3 — what leakage actually is.** The leakage of one layer is exactly
a weighted variance: take the quantity g(y) = Σ_d f_d(β)·n(y;d,c') and ask how
much it varies *within* each cost class. Uniform within classes → zero leakage
→ proxy exact. This upgrades the original paper's empirical "the n(x;d,c)
don't vary much" plots into an identity, with one twist that matters: what
must concentrate is a **β-weighted combination** of the n's, not each entry
separately. That twist is our best explanation for an old mystery (see the
fitted-shape paradox below). Verified to 4.5e-15.
→ [exp 008](experiments/008_sampled-leakage-predictor/README.md)

## What the experiments established, as a story

**1. Is the game worth playing? Yes.** On every graph family we tested (dense
and sparse Erdős–Rényi, Barabási–Albert, Watts–Strogatz, 3-regular — 420
instances), proxy-set parameters beat the random-balanced-partition baseline
on *every single instance*, and the compression itself is nearly lossless at
p=1: the exact-compression proxy lands within ~0.03 approximation ratio of
the best any parameters could do.
→ [exp 002](experiments/002_baselines-and-headroom/README.md)

**2. What makes compression fail? Density — not "how ER-like" the graph is.**
The surprise of Phase 1: family-average leakage tracks *edge count* with
correlation 0.97. Dense ER(0.5) is the *worst*-compressing family we tested;
3-regular the best. The community's intuition "the proxy is an ER trick"
points at the wrong variable. → [exp 003](experiments/003_leakage-vs-overlap/README.md),
sharpened by the small-angle law λ ≈ const·β·γ²·m
→ [exp 007](experiments/007_leakage-anatomy/README.md)

**3. Does leakage predict *parameter-setting* quality, not just state
fidelity? Yes at p=1.** Families ranked by leakage are ranked by regret
(how far below the grid ceiling the proxy's chosen parameters land), Spearman
ρ = 0.86–0.96. This was the go/no-go gate for the whole theory-led framing —
it passed. At p=3 along proxy-chosen ramps, both leakage and regret flatten
across families, so there's no ranking signal there (not a contradiction —
there's just nothing to rank). → [exp 004](experiments/004_gate-leakage-vs-regret/README.md)

**4. The analytical formula is strangely robust — and its one failure is
diagnosable.** The Sud et al. analytical N (with an effective edge
probability) sets parameters *better* than exact compression on sparse and
non-ER families — regret 0.01–0.02 — because the smooth formula acts like a
regularizer. Its predicted *values* are garbage (badly inflated norms), but
its *argmax location* — the thing parameter setting actually uses — is
excellent everywhere except dense ER(0.5), where a spurious large-β peak
appears. We tried two "filter" fixes (norm sanity checks, physicality caps);
both fail as universal recipes because the model's calibration breaks long
before its argmax moves. Verdict: don't trust the analytical proxy's numbers,
trust its peak — except on dense graphs.
→ [exp 005](experiments/005_norm-filtered-paper-proxy/README.md),
[exp 006](experiments/006_physicality-filter/README.md)

**5. Compression error and model error really are independent axes.** The
spurious dense-ER peak (a model-error artifact) is invisible to leakage (a
compression-error measure) — it sits at the 38th percentile of leakage.
You need both diagnostics; neither substitutes for the other.
→ [exp 007](experiments/007_leakage-anatomy/README.md)

**6. You don't need the whole graph — 5 samples per cost class suffice.**
The leakage predictor estimated from S=5 stratified samples matches the exact
value to ~3% median error, uniformly across families and angles. This is what
makes everything above *practical* at scales where you can't enumerate 2^n
bitstrings. → [exp 008](experiments/008_sampled-leakage-predictor/README.md)

**7. Why does the proxy degrade at depth? Wrong aim, not wrong size.** The
true QAOA trajectory under ramp schedules is astonishingly low-dimensional
(2–4 dimensions capture 99% of it — even at p=20). So the cost-class subspace
(m+1 dimensions!) is more than big enough; when the proxy degrades, it's
because the *fixed* cost-class frame points slightly away from where the
trajectory actually lives. → [exp 009](experiments/009_trajectory-pca/README.md)

**8. Everything survives scale-up.** At n=16–18 (where 2^n is a quarter
million states): regret stays flat in n, value-added stays positive, and —
the practical headline — the *sampled*-N proxy (S=10) matches exact-N
parameter choices to ~0.01 AR, sometimes better. The practical recipe is:
sampled N + empirical P, no filters needed.
→ [exp 010](experiments/010_scaleup-ranking/README.md)

**9. Depth is predictable too.** Accumulated leakage grows *linearly* in p for
fixed ramps — measured ratio 1.505 vs the 1.5 Theorem 3 predicts — and only
mildly with n. Even at p=30, n=20, small ramps keep 0.71–0.81 overlap with
the true state. And notably: the proxy's parameter-setting success at depth
operates at fidelities well *below* 1 — what survives depth is the argmax,
not the state. → [exp 011](experiments/011_depth-scaling/README.md)

**10. What does it all cost? (New, 2026-07-03.)** Honest wall-clock on an
A100: the pipeline's bottleneck is *obtaining N*, not running the proxy.
Exact empirical N is O(4^n) — 21.8 s at n=20, hopeless past n≈24. Sampled N
is 50× cheaper. The proxy sweep itself is nearly independent of n (it scales
with edge count squared). The honest caveat for the paper: at sizes we can
simulate, GPU statevector brute force is just as fast as the proxy — the
proxy's real value is at scales beyond simulation, or when each real-QAOA
evaluation is expensive (actual hardware).
→ [exp 013](experiments/013_timing-benchmark/README.md)

## The open mystery we're attacking next: the fitted-shape paradox

The G-RIPS triangle and normal proxies were *fitted* to match the empirical
N(c';d,c) — and fits with *lower* mean-squared error set parameters *worse*.
That's backwards, unless the error that matters isn't entrywise. Theorem 3
says exactly that: the proxy update only sees N through the β-weighted
combination Σ_d f_d(β)·N(c';d,c), where cancellations and weightings make
some entries matter enormously and others not at all. Experiment 014 (now
running as Slurm job 11575420) tests this directly: perturb N in ways that
are large entrywise but small in the weighted norm (and vice versa), and
refit the triangle/normal shapes, checking whether the weighted model error —
not the entrywise MSE — predicts parameter-setting regret. If it does, the
paradox dissolves and the paper gains a clean prescription: *fit proxies in
the norm the algorithm actually uses.*

**A discovery from E014's shakedown run, in plain language.** The proxy
reports a predicted average cut value for each candidate (γ, β), and we pick
the candidate with the biggest prediction. It turns out that if the N you
feed the proxy is even slightly wrong — 1% off, even in a direction the
dynamics provably can't see — the predicted values at *some* corner of the
parameter grid explode (we saw a prediction of ~500,000 on a graph with 21
edges, whose best cut is therefore at most 21). The reason: with a wrong N,
the compressed state's total probability can grow past 1, and the "expected
cost" of an inflated state looks huge. With the exact N this can never happen
— Theorem 1 guarantees the total probability never grows — and we watched
that hold: exact N kept total weight at 0.978 while every wrong N inflated
it, up to tens of thousands. The candidate fix is embarrassingly simple:
divide the prediction by the total weight, i.e. ask "what is the average cost
of the state we actually have?" On the shakedown instance that one division
made the theory-invisible perturbations exactly harmless (zero regret, as
Theorem 3 says) and made regret track the weighted error cleanly. One
wrinkle keeps us honest: the paper's own analytical proxy got slightly
*worse* after the division on that instance — so no victory declaration until
the 150-instance run reports. Earlier experiments (005/006) tried *rejecting*
suspiciously inflated predictions and failed; *dividing* by the inflation may
be the recipe that works.

## Where the paper stands

Full IEEE-format draft exists (`papers/OverleafPaper/qce2027_paper.tex`,
branch ClaudeResearch) with the theory, experiments 001–011, and the E013
pipeline-cost subsection; the 4-lens review fixes are applied and pushed.
Remaining before the August arXiv target: resolve E014 and revise §6
accordingly, fold in the systematic literature pass (scan running, report →
`lit_scan_2026-07-03.md`), figures/tables polish, author list.
