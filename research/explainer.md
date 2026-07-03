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

## The fitted-shape paradox: solved (exp 014, 2026-07-03)

The G-RIPS triangle and normal proxies were *fitted* to match the empirical
N(c';d,c) — and fits with *lower* mean-squared error set parameters *worse*.
That's backwards, unless the error that matters isn't entrywise. Theorem 3
says exactly that: the proxy update only sees N through the β-weighted
combination Σ_d f_d(β)·N(c';d,c), where cancellations and weightings make
some entries matter enormously and others not at all.

**11. The paradox dissolves under the right error metric (exp 014,
150 instances, 5 families, n=12–14).** We corrupted the exact N with
same-sized errors pointing in different "directions": one direction the
dynamics amplifies, one it provably cannot see, one in between, one random.
Entrywise, these corruptions are identical in size — up to 50% of N. Result:
the invisible corruptions cost *exactly nothing* (regret equal to the exact-N
baseline, even at 50% error), while the amplified direction at just 1% error
already costs 0.13 of approximation ratio. Across everything we ran, the
weighted error predicts regret at rank-correlation 0.74; entrywise MSE
manages 0.38 — and on the actual triangle-vs-normal comparison MSE gets the
order *backwards* in 105 of 150 instances (the normal fit has lower MSE but
higher regret). The old G-RIPS observation was real, and it was never a
paradox — we were measuring fit quality with the wrong ruler.
→ [exp 014](experiments/014_fitted-shape-paradox/README.md)

**12. But the tempting fix fails — the shapes themselves are the problem.**
If MSE is the wrong ruler, surely refitting the triangle/normal shapes with
the *right* ruler (the weighted error) fixes parameter setting? No: the
refit barely moves the triangle (better in 77 of 150, worse in 67) and hurts
the normal. The reason is honest and final: even the *best possible*
triangle or normal keeps ≥ 92% relative error in the weighted norm — these
two-to-four-parameter shapes simply cannot imitate how the real N drives the
dynamics. Shape fitting didn't fail because of a bad objective; it failed
because the shape families are too poor. This strengthens the paper's actual
recommendation: estimate N by *sampling bitstrings* (a few per cost class,
exp 010) — sampling has no shape assumption to get wrong.

**13. A subtlety about "inflated" predictions that sharpens experiments
004–006.** With any slightly-wrong N, the proxy's predicted values explode
somewhere on the parameter grid (predictions of ~500,000 on 21-edge graphs
whose best cut is 21): wrong N lets the compressed state's total probability
grow past 1, which Theorem 1 forbids for the exact N — and indeed exact N
stayed ≤ 0.984 on all 150 instances while every wrong N inflated, sometimes
by 10⁵. The obvious cure — divide each prediction by its state's total
weight, "average cost of the state you actually have" — is the right
*measuring instrument* (results 11–12 use it), but as a parameter-setting
recipe it *backfires*: it slightly hurts even the exact N, and it badly
hurts the paper's analytical proxy (regret 0.05 → 0.16). Translation: the
analytical proxy's raw numbers are meaningless as *values* (that we knew),
but *where they are large* is genuinely informative — its inflation tends to
sit on top of good parameters (except the known dense-graph artifact), so
dividing the inflation out throws away signal. Practical bottom line is
unchanged: sampled N + empirical cost distribution + the ordinary objective.

## Where the paper stands

Full IEEE-format draft exists (`papers/OverleafPaper/qce2027_paper.tex`,
branch ClaudeResearch) with the theory, experiments 001–011, and the E013
pipeline-cost subsection; the 4-lens review fixes are applied and pushed.
E014 is resolved (results 11–13 above) and now needs to be *written into*
§6, including a likely figure (regret vs corruption size for the four
directions — the visual proof that error direction, not size, is what
matters). Also remaining before the August arXiv target: fold in the
systematic literature pass (scan running, report →
`lit_scan_2026-07-03.md`), figures/tables polish, author list.
