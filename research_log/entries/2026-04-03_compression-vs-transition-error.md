# Compression vs Transition Error: Why Non-ER Fails and Why MSE ≠ Performance

**Date:** 2026-04-03
**Paper section:** discussion
**Tags:** theory, homogeneity, non-ER, proxy-failure, homodist-MSE
**Status:** analysis (no new experiments)

## Motivation

Two major negative results need deeper explanation:
1. Non-ER is closed as a proxy research path — SampN+EmpP fails on all BA/WS
2. Homodist MSE does not predict QAOA performance — fitted proxies with low MSE
   produce bad parameters

The existing explanations ("homogeneity breaks down", "errors amplify") are
descriptive but don't identify the fundamental mechanism.

## Two Independent Sources of Proxy Error

The proxy compresses a 2^n-dimensional state vector into an (m+1)-dimensional
vector (one amplitude per cost class). This introduces two distinct error types:

### Error (1): Compression Error

Information lost by projecting 2^n amplitudes onto m+1 cost classes. The true
QAOA state has *different* amplitudes for bitstrings with the same cost:

    ψ_true(x) ≠ ψ_true(y)  even when  c(x) = c(y)

The proxy forces ψ_proxy(x) = Q(c(x)) for all x — a rank-(m+1) approximation
of the full state. The compression error is the norm of the component
orthogonal to the homogeneous subspace:

    ε_compress = ||ψ_true - Π_homo(ψ_true)||

where Π_homo projects onto the subspace of cost-class-uniform states.

This error depends on the *graph structure*, not on the proxy's N(c';d,c).
No choice of N can reduce it.

### Error (2): Transition Error

Error in the transition matrix within the compressed subspace. The proxy uses
N(c';d,c) to propagate amplitudes between cost classes. If N is wrong, the
dynamics within the homogeneous subspace are wrong:

    ε_transition = ||M_proxy - M_true_projected||

where M_true_projected is the transition matrix the true dynamics would induce
on the homogeneous subspace. This is what homodist MSE measures.

## Why Non-ER Fails (Compression Error Dominates)

For ER graphs, the random structure means vertices are statistically
interchangeable. Two bitstrings with the same cut value involve cutting
statistically similar edges (no hubs, no clusters). Their neighborhoods
under the mixer are similar, so their amplitudes stay close after each
layer. Compression error is small.

For BA graphs, hubs create heterogeneity: a bitstring that cuts edges around
a high-degree hub behaves very differently under the mixer than one that cuts
peripheral edges, even if both have the same total cost. For WS graphs,
cluster structure creates similar heterogeneity.

**This is why SampN+EmpP fails on non-ER despite having near-zero transition
error.** The empirical homodist is the *correct* average N(c';d,c), so
transition error is minimal. But the compression error is large because
the true state rapidly leaves the homogeneous subspace. No improvement to
N can compensate — the proxy formalism itself is insufficient for these graphs.

**Why the gap worsens with depth:** Each layer introduces new within-cost-class
variance. Even if the initial state is perfectly homogeneous (uniform
superposition), after one layer of phase + mixer, bitstrings in the same cost
class acquire different amplitudes based on their local graph structure. After
p layers, the accumulated compression error dominates.

**Why the gap worsens with n (on non-ER):** Larger BA/WS graphs have more
pronounced structural heterogeneity (wider degree distributions, more distinct
community structure), increasing per-layer compression error.

## Why Homodist MSE ≠ QAOA Performance (Wrong Loss Function)

Even on ER where compression error is small, fitted proxies (Triangle, Gaussian)
fail despite low homodist MSE. The reason is that homodist MSE is the wrong
optimization target.

The proxy objective landscape is:

    f(γ,β) = Σ_{c'} P(c') |[M_p · ... · M_1 · Q_0]_{c'}|² · c'

where each M_l depends on (γ_l, β_l) and N(c';d,c). The landscape topology
(location and relative heights of optima) determines which parameters the
optimizer finds.

Minimizing ||N_proxy - N_true||_F (homodist MSE) minimizes the per-element error
in each M_l. But the landscape depends on the *product* M_p · ... · M_1, and
small per-element errors can steer this product in arbitrary directions.

More precisely, the error in the landscape gradient is:

    ∇f_proxy - ∇f_true ≈ Σ_l [∂f/∂M_l] · δM_l

where δM_l comes from the homodist error. The perturbation that minimizes
||δM_l||_F (MSE) is generically NOT the perturbation that minimizes
||∇f_proxy - ∇f_true|| (landscape gradient error). These are different
optimization problems with different solutions.

**Why PaperProxy succeeds despite being "wrong" for non-ER**: Its errors are
*structurally correlated* between N and P, because both are derived from the
same probabilistic model (binomial/multinomial over independent edges). The
error in N and the weighting by P partially cancel in the objective function.
Fitted proxies have errors that are uncorrelated with P, so they don't cancel.

**Why smoothing destroys performance**: Smoothing reduces high-frequency
transition error but (a) introduces low-frequency transition error in its
place, (b) doesn't touch compression error, and (c) destroys the fine
structure in N that determines the landscape's critical points.

## Testable Predictions

1. **Within-cost-class amplitude variance**: Compute Var[ψ_true(x) | c(x)=c']
   after each QAOA layer. Predict: higher for non-ER than ER, grows with depth
   for both, rate of growth correlates with proxy failure severity.

2. **Per-instance proxy performance vs homogeneity**: Among ER instances, those
   with higher within-cost-class variance (less homogeneous) should have worse
   proxy performance. This would confirm that compression error, not transition
   error, is the binding constraint.

3. **Proxy with variance tracking**: A "second-order proxy" that tracks both
   the mean and variance of amplitudes within each cost class could potentially
   correct for compression error. Each layer would propagate both Q(c') and
   σ²(c'), using the variance of n(x;d,c) within each cost class.

## Implications for the Paper

This framework provides a clean theoretical contribution:

1. **Characterize the proxy as a dimensionality-reduction method** with two
   error sources (compression + transition), not just a "heuristic."

2. **Explain when the proxy works**: when compression error is small (ER) and
   transition error is small (analytical or empirical N with consistent P).

3. **Explain the non-ER negative result constructively**: it's not that our
   implementation is bad — it's that cost-class compression is provably lossy
   for structured graphs. This scopes the proxy method's applicability.

4. **Frame homodist MSE ≠ performance as a misaligned loss function problem**:
   the right metric would be landscape fidelity, not distribution fidelity.

## Low-Rank Framing (Key Insight)

The compression error framework maps directly onto low-rank matrix approximation.
View n(x;d,c) as a matrix with rows indexed by x (2^n bitstrings) and columns
indexed by (d,c) pairs. The homogeneous approximation replaces this with
N(c(x);d,c), which has at most m+1 distinct rows — i.e., it is a rank-(m+1)
matrix whose row space is spanned by cost-class indicator vectors.

By Eckart-Young, the best rank-k approximation is given by the truncated SVD.
The proxy quality therefore depends on:
1. **Singular value decay**: If the top m+1 singular values capture most of the
   Frobenius norm, the rank-(m+1) approximation is good (proxy works).
2. **Basis alignment**: Whether the top singular vectors align with cost-class
   indicators. Even if rank-(m+1) suffices, the proxy uses a specific basis
   (cost classes), not the optimal SVD basis.

**Predictions:**
- ER graphs: fast singular value decay, good alignment → proxy works
- BA/WS graphs: slow decay, poor alignment → proxy fails
- The singular value gap at rank m+1 is a quantitative, a priori predictor of
  proxy quality for any graph family

**Constructive direction:** If the SVD basis doesn't align with cost classes,
a higher-rank proxy using SVD basis vectors (instead of cost-class indicators)
could extend the method to structured graphs. This would be a novel contribution
beyond the scope of Sud et al.

## Other Speculative Directions

- **Second-order proxy**: Track variance alongside mean amplitude per cost class.
  Propagate both through each layer using the variance of n(x;d,c).

- **Selective decompression**: Use the proxy for most cost classes but track
  individual amplitudes for the few cost classes with highest within-class
  variance.

---

*Analysis entry, 2026-04-03*
