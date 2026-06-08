# Julia vs QOKit statevector agreement

**Date:** 2026-06-08
**Script:** `scripts/analysis/statevector_error_analysis.jl`
**Context:** Quantifies the numerical agreement behind the CI test
`test/test_qaoa_simulation.jl` ("qaoa_statevector matches QOKit"), which
failed because its tolerance (`atol=1e-18`) was set far below double-precision
rounding. This report establishes what the real error magnitude is, so the
tolerance can be chosen on evidence.

## Setup

For each random instance we build a G(n, 0.5) Erdős–Rényi MaxCut graph, draw
random QAOA parameters γ, β ∈ [0, 2π), and compare:

- `jl_state` — `JuliaQAOA.qaoa_statevector` (Julia FUR simulation), against
- `py_state` — QOKit's reference Python FUR backend
  (`grips.QAOA_simulator.get_state`, `simulator_name="python"`).

Both produce a length-`2^n` complex statevector. We measure, per amplitude,
the absolute error `|Q_jl − Q_py|` and relative error `|Q_jl − Q_py| / |Q_py|`,
and per state the 2-norm difference `‖Q_jl − Q_py‖₂`. Sweep:

| n  | p | instances |
|----|---|-----------|
| 8  | 3 | 50        |
| 10 | 5 | 20        |
| 12 | 8 | 10        |

(The CI test itself uses n=8, p=3, seeds 0–4.)

## Results

**Per-amplitude absolute error** `|Q_jl − Q_py|`:

| statistic | value      |
|-----------|------------|
| mean      | 1.693e-17  |
| median    | 1.362e-17  |
| max       | 2.779e-16  |

**Per-amplitude relative error** `|Q_jl − Q_py| / |Q_py|`:

| statistic | value      |
|-----------|------------|
| mean      | 6.998e-16  |
| median    | 3.487e-16  |
| max       | 9.830e-13  |

**Whole-state 2-norm difference** `‖Q_jl − Q_py‖₂`:

| statistic | value      |
|-----------|------------|
| mean      | 4.770e-16  |
| median    | 3.522e-16  |
| max       | 1.642e-15  |

**Per-configuration worst case:**

| config           | max abs err | max ‖Δ‖₂  |
|------------------|-------------|-----------|
| n=8,  p=3 (50)   | 2.779e-16   | 1.642e-15 |
| n=10, p=5 (20)   | 9.813e-17   | 1.137e-15 |
| n=12, p=8 (10)   | 6.069e-17   | 1.380e-15 |

## Interpretation

- Absolute amplitude errors sit at the level of a single double-precision ULP
  for amplitudes of magnitude ~0.01–0.1 (`eps(Float64) ≈ 2.2e-16`). This is the
  expected consequence of the two backends accumulating the same arithmetic in
  a different order (different loop structure / BLAS), not a correctness
  difference.
- The whole-state 2-norm difference never exceeds **1.64e-15** across all 80
  instances, including deeper circuits (p=8) and larger states (n=12). It does
  **not** grow with n or p in this range.
- The large *relative* max (9.8e-13) is an artifact of near-zero amplitudes:
  when `|Q_py|` is itself ~1e-16, the ratio inflates even though the absolute
  error is a single ULP. The absolute and 2-norm measures are the meaningful
  ones here.

## Tolerance choice

The test's `@test jl_state ≈ py_state atol=1e-12` is justified: the observed
worst-case `‖Δ‖₂` (1.64e-15) clears it by ~3 orders of magnitude, so the test
is robust to platform/BLAS variation, while still being ~3 orders of magnitude
tighter than the next looser comparison in the file. The original `atol=1e-18`
was ~1000× below a single ULP and could not pass on any platform.
