# 021 — The section-4 theory extends verbatim to integer-weighted MaxCut

**Question:** The paper's Limits scope everything to unweighted MaxCut. Does
the theory chain (neighbor-sum lemma, variance reduction, codegree form, the
four moment identities, the quadratic conditioning bound, the cubic leakage
law) survive edge weights?

**Answer: Yes, verbatim, for integer weights, with the natural weighted
replacements; and the unweighted restriction was never about the theory:
continuous weights dissolve the compression itself by collapsing every cost
class to a bitstring-complement pair.**

## Replacements (all verified to 1e-9 on 140 instances, 7 families x
n in {12,14} x 10, iid weights in {1,2})

- m -> W (total weight); Var(S) = m2 = sum w^2; q4 = sum w^4
- codegree A_jk -> A^w_jk = sum_i w_ij w_ik
- triangles tau -> tau_w = sum_triangles w_ab w_bc w_ca
- 4-cycles c4 -> c4_w = sum_4cycles w_ab w_bc w_cd w_da

Then, exactly as in the unweighted case:
- (L1w) sum_i c(x xor e_i) = (n-4) c(x) + 2W on every bitstring
- (L3w) E[T] = 2 m2, Var(T) = 2 sum_{j!=k} (A^w_jk)^2
- E[S^3] = 6 tau_w; Var(S^2) = 2 m2^2 - 2 q4 + 24 c4_w
- Cov(T,S) = 6 tau_w; **Cov(T,S^2) = Var(T)** (the operative identity again)
- law of total variance closes: Var(T) = V2 + Var(E[T|c])

## Measured

- **Quadratic conditioning bound capture: mean 0.975, range 0.929-0.996**
  over the 140 weighted instances (zero violations) --- the same 97% as the
  unweighted E019.
- **Cubic law at small angles: lambda_1 / [(|beta| gamma^2 / 8) sqrt(V2)] =
  0.998 +- 0.001** (range 0.996-1.001) at (gamma, beta) = (0.05, 0.05).

## The structural boundary

With iid uniform(0,1) weights (n=12, one instance per family): the number of
distinct cost values is exactly 50% of 2^n and the largest cost class has
size 2 on every family --- each bitstring shares its cost only with its
complement. The homogeneous subspace is then (essentially) the whole Hilbert
space: the compression compresses nothing, its "exactness" is vacuous, and
the proxy's O(m^2 n) evaluation advantage disappears (the number of classes
is exponential). So the proxy requires cost degeneracy; integer (or small
rational) weights preserve it, generic continuous weights destroy it. This
is a boundary of the OBJECT, not of the theory.

## Paper impact

- New remark in section 4: the whole section extends to integer-weighted
  MaxCut with the replacements above (verified, E021); continuous weights
  dissolve the compression itself.
- Limits updated: "unweighted MaxCut" scope now applies to the experiments;
  the theory is stated weighted-ready.

Reproduce: `JULIA_NUM_THREADS=auto julia --project research/experiments/021_weighted-maxcut/run.jl` (~4 min).
