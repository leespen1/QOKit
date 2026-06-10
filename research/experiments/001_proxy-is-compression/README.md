# 001 — Is the homogeneous proxy exactly a subspace compression?

**Question:** Is the homogeneous proxy iteration (run with the same-instance
empirical N(c'; d, c)) numerically identical to evolving the true statevector
one QAOA layer at a time and orthogonally projecting onto cost classes after
each layer?

**Answer: Yes — identical to machine precision (max deviation 2.5×10⁻¹⁶,
threshold 10⁻¹⁰), so "proxy = compression" (Theorem 1) is exact, not an
approximation.**

## Method

For 15 Erdős–Rényi(0.5) MaxCut instances (5 each at n = 8, 10, 12; p = 3 with
small/moderate/large angles), we computed the proxy trajectory two completely
independent ways:

1. **Distribution route:** build the empirical homogeneous distribution
   N(c'; d, c) with `get_homogeneous_distribution_from_costs_direct`, then run
   the proxy recursion (`QAOA_proxy_basic`, and `QAOA_proxy_single` as a second
   implementation).
2. **Projection route:** evolve the full 2^n statevector one layer at a time
   (`apply_phase_gate!` + `apply_x_mixer!`) and after each layer replace every
   amplitude by the mean over its cost class
   (`compressed_qaoa_trajectory` in `src/subspace_compression.jl`).

If the proxy is exactly the orthogonal compression of QAOA onto the
cost-class subspace, the two must agree on every attained cost class at every
layer. We also confirmed the `pi_units=true` convention (used by the figure
scripts) reproduces the raw-radians result.

## Result

Max absolute deviation over all instances, layers, and classes:

| comparison                         | max \|Δ\|   |
|------------------------------------|-------------|
| projection vs `QAOA_proxy_basic`   | 1.3×10⁻¹⁶  |
| projection vs `QAOA_proxy_single`  | 1.1×10⁻¹⁶  |
| radians vs `pi_units=true`         | 2.5×10⁻¹⁶  |

## Caveats

- **Unattained cost classes:** the proxy seeds *all* classes 0..m with the
  uniform amplitude 1/√2ⁿ at layer 0, including costs no bitstring attains;
  the projection gives 0 there. These entries are inert (the empirical N has
  zero rows/columns for them), so both evolutions hold them at exactly 0 from
  layer 1 onward. Agreement is therefore checked on attained classes at layer
  0 and on all classes afterward.
- This exactness holds only for the **same-instance empirical N**. The
  analytical PaperProxy N (and any fitted N) is *not* the compression of any
  specific instance — its deviation is a separate "model error" on top of the
  compression error. Separating those two is the subject of Phase 1.

## Reproduce

```
julia --project research/experiments/001_proxy-is-compression/run.jl
```

Seed 20260610 fixed in the script; runs in under a minute on a laptop.
Unit-test version (n = 6, 8): `test/test_subspace_compression.jl`.
