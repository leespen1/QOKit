# Theory notes: where Theorems 1–3 come from, and prior art

*Written 2026-06-11 in response to Spencer's question about the theoretical
foundations. This is the basis for the paper's §3 and Related Work.*

## The one reframing everything follows from

Sud et al. (arXiv:2211.09270, PRResearch 6, 023171) derive the proxy by a
*substitution*: assume amplitudes equal within cost classes, propagate one
exact layer (sum-of-paths), then "replace n(x;d,c) with N(c';d,c)." They state
explicitly that after this replacement "this evolution is no longer restricted
to unitary evolution" and Q is only an "analogue of amplitude." Validity is
checked empirically (overlap plots, stddev/mean heatmaps).

Our reframing asks: *what linear map on the full 2^n-dimensional space does the
recursion implement when N is the same-instance empirical class average?* The
answer is exact, not approximate, and everything else is elementary linear
algebra from there.

## Theorem 1 (proxy = orthogonal compression) — one line of algebra

With N(c';d,c) = (1/M_{c'}) Σ_{x∈S_{c'}} n(x;d,c) (what
`get_homogeneous_distribution_from_costs` computes), plug into the proxy update
and swap summation order:

    Q_ℓ(c') = Σ_{d,c} f_d(β) e^{-iγc/2} Q_{ℓ-1}(c) · (1/M_{c'}) Σ_{x∈S_{c'}} n(x;d,c)
            = (1/M_{c'}) Σ_{x∈S_{c'}} [ Σ_{d,c} f_d(β) e^{-iγc/2} Q_{ℓ-1}(c) n(x;d,c) ]
            = (1/M_{c'}) Σ_{x∈S_{c'}} q_ℓ(x),

where q_ℓ(x) is the *exactly evolved* amplitude (Sud et al.'s own Eq. before
substitution). So one proxy step = evolve exactly, then replace each amplitude
by its cost-class mean. "Replace by class means" is precisely the orthogonal
projection P onto span{|c⟩ = M_c^{-1/2} Σ_{x∈S_c}|x⟩} (the classes partition
{0,1}^n, so the |c⟩ are orthonormal). Hence the proxy trajectory is
(PU_p)···(PU_1)|+⟩^n in suitable coordinates: an exact compression, verified
to machine precision in experiment 001. The non-unitarity Sud et al. flag is
not a defect — it is exactly the contractivity of a projection, and the norm
loss per layer is a meaningful, measurable error (the leakage).

## Theorem 2 (telescoping bound) — standard perturbation argument

ψ_ℓ = U_ℓ ψ_{ℓ-1} and φ_ℓ = P U_ℓ φ_{ℓ-1}. Then

    ψ_ℓ − φ_ℓ = U_ℓ(ψ_{ℓ-1} − φ_{ℓ-1}) + (I−P) U_ℓ φ_{ℓ-1},

take norms; unitarity of U_ℓ preserves the first term, the second is λ_ℓ.
Induction gives ‖ψ_p − φ_p‖ ≤ Σ_ℓ λ_ℓ. This proof pattern is textbook
(Duhamel/telescoping; the same skeleton as Trotter error bounds and a
posteriori bounds in projection-based model order reduction). Nothing novel in
the technique — the contribution is *applying* it with a λ_ℓ that costs one
statevector layer + an O(2^n) projection to measure, and showing empirically
(exp 003) that it is ~4× from tight and that Σλ is a near-functional fidelity
predictor.

## Theorem 3 (leakage = weighted within-class variance) — law of total variance

UP|c'⟩ has computational-basis amplitudes proportional to
g_{c'}(y) = Σ_d f_d(β) n(y;d,c'). (I−P) subtracts class means, so

    λ(c')² ∝ Σ_c Σ_{y∈S_c} |g_{c'}(y) − mean_{S_c}(g_{c'})|²
           = Σ_c M_c · Var_{y∈S_c}[g_{c'}],

i.e. exactly the class-size-weighted within-class variance of g. This upgrades
Sud et al.'s empirical §5 check ("stddev of n(x;d,c) around its class mean is
small for dominant terms") into an identity, with one refinement that matters:
the relevant concentration is of the f_d(β)-weighted combination, not
entrywise — which is the candidate explanation for why entrywise-MSE proxy
fitting failed (hypothesis H-fit) and for the small-angle regime.

## Prior art map (for Related Work)

- **Hogg (2000), quant-ph/0006090 etc.**: the original "mean-field"/homogeneous
  model the proxy generalizes; ansatz-level, no projection framing or bounds.
- **Sud, Hadfield, Rieffel, Tubman, Hogg (2022/2024)**: the proxy + heuristic;
  substitution derivation; explicitly non-unitary "analogue amplitudes";
  validation entirely empirical. No P, no λ, no bounds.
- **Shaydulin, Hadfield, Hogg, Safro, "Classical symmetries and the QAOA"
  (arXiv:2012.04713)**: *exact* invariant subspaces from cost-function
  symmetries compatible with the mixer; Sud et al. themselves note Perfect
  Homogeneity is stronger than symmetry. Our work is the approximate-invariance
  generalization of this line: the cost-class subspace is generally NOT a
  symmetry subspace, so invariance fails by a measurable amount (leakage).
- **"Symmetries and Dimension Reduction in QAOA" (arXiv:2309.13787)** and
  Dicke-subspace QAOA work: same exact-symmetry theme.
- **Lumping / coarse-graining / projection-based model order reduction**
  (Markov-chain lumpability; e.g. arXiv:1104.1025 on hidden-Markov
  coarse-graining error, arXiv:2512.11974 on lumping limits; standard MOR
  a posteriori bounds): the mathematical machinery for "project linear
  dynamics onto a partition subspace and bound the error by accumulated
  residual" is classical. Cost classes = a lumping partition of {0,1}^n;
  our Theorems 2–3 are this machinery instantiated for unitary QAOA layers.
- **DiezValle et al. pseudo-Boltzmann (cited by Sud)**: p=1 amplitude-as-
  function-of-cost structure; orthogonal motivation.

**Claimed-new in our paper (pending final lit pass):** (a) the identification
of the empirical-N proxy as *exact* orthogonal compression (Thm 1) rather than
a heuristic substitution; (b) per-layer leakage as a cheap measured diagnostic
with the telescoping bound for QAOA (Thm 2); (c) the variance identity with
f_d(β) weighting (Thm 3); (d) the empirical program: leakage ⇔ fidelity ⇔
regret rankings, the density (not ER-ness) finding, model-error/compression-
error separation, argmax-robustness of the analytical N, and the artifact
anatomy. The *proof techniques* are standard and will be presented as such.
