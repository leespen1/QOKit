# 013 — Timing: what does each stage of proxy parameter setting cost?

**Question:** What is the honest wall-clock cost of every stage of the proxy
parameter-setting pipeline (obtaining N, sweeping schedules, real-QAOA
reference evaluations), and how does each scale with n?

**Answer: The pipeline is dominated by obtaining N, not by running the proxy —
exact N costs O(4^n) (21.8 s at n=20 on an A100, prohibitive beyond n≈24)
while sampled N (S=10) costs 0.44 s — and at classically simulable sizes the
proxy buys no wall-clock advantage over GPU statevector grid search (p=1 sweep:
0.38 s proxy vs 0.48 s brute-force ceiling at n=20; even at p=20 ramps the
comparison stays even, est. 16 s brute force vs 13.4 s sweep + 0.4 s sampled
N); the proxy's value is the polynomial-vs-2^n asymptotic, the reuse of one N
across every depth and schedule grid, and settings where each real-QAOA
evaluation is expensive, not speed at n≤20.**

## Method

One G(n, 0.5) instance per n ∈ {12, 14, 16, 18, 20} (seed 20260612 + n), on
one HPCC node (NVIDIA A100-SXM4-80GB; AMD EPYC 7713, 8 Julia threads). Each
operation is timed as the median of 3 reps after a warmup call (JIT excluded;
300 s budget per op). Stages:

- **exact_N** — empirical N(c';d,c) via the GPU direct homodist kernel, O(4^n).
- **sampled_N** — stratified sampled N, S=10 per cost class (CPU, threaded).
- **proxy_sweep_p1/p3/p20** — the full parameter-setting decision (sweep +
  expectation + argmax) via `QAOA_proxy_multi`: 40×40 (γ,β) grid at p=1;
  8⁴ = 4096 linear-ramp endpoint schedules at p=3 and p=20. CPU BLAS —
  the proxy never needs the GPU.
- **sv_single_p1/p20** — one real-QAOA expectation (GPU batched statevector).
- **ceiling_p1/p3** — brute-force grid ceiling: 1600 (p=1) / 4096 (p=3)
  statevector evaluations.

Grids and bounds identical to experiments 002/004/010. GPU timings exclude
the one-time cost upload (steady-state regime, matching how the ceilings were
computed in 002/004/010).

## Result

Wall-clock seconds (n=20, m=103, A100):

| stage | seconds |
|---|---|
| exact_N (GPU, O(4^n)) | 21.8 |
| sampled_N (S=10, CPU) | 0.44 |
| proxy_sweep_p1 (1600 schedules) | 0.38 |
| proxy_sweep_p3 (4096 ramps) | 2.3 |
| proxy_sweep_p20 (4096 ramps) | 13.4 |
| sv_single_p1 | 0.00034 |
| sv_single_p20 | 0.0040 |
| ceiling_p1 (1600 statevector evals) | 0.48 |
| ceiling_p3 (4096 evals) | 2.7 |

Full table for all n in `results.csv`. Three observations:

1. **Obtaining N is the bottleneck, and sampling fixes it.** exact_N grows
   ~13–17× per +2 qubits (0.094 s → 1.25 s → 21.8 s for n=16→18→20),
   consistent with O(4^n): extrapolating, ~6 min at n=22, ~1.5 h at n=24.
   sampled_N grows ~4× per +2 qubits (O(2^n) to touch the cost array) with a
   tiny constant — 50× cheaper than exact at n=20, and E3.1 (exp 010) showed
   S=10 sampled N matches exact-N parameter choices within ~0.01 AR.
2. **The proxy sweep itself is cheap and nearly n-independent.** Its cost
   scales with m² (BLAS mat-mat on (m+1)-dim states), not 2^n: the p=1 sweep
   is 0.037 s at m=36 and 0.38 s at m=103. Depth enters linearly
   (p=3: 2.3 s → p=20: 13.4 s for the same 4096 schedules).
3. **At simulable sizes, brute force is as fast as the proxy.** A single GPU
   statevector expectation is kernel-launch-bound at these sizes (~0.2–0.3 ms,
   nearly flat in n up to 20), so the 1600-point p=1 ceiling (0.48 s) costs
   about the same as the proxy sweep (0.38 s) *before* counting the cost of
   obtaining N. The pipeline's cost advantage opens only as 2^n outruns launch
   overhead, around n ≈ 22–24 (extrapolated) — and there, N must come from
   the analytical formula or sampling, since exact N is O(4^n). This is the
   honest framing for the paper's practical-consequences section.

## Caveats

- One instance per n; timings are stable across reps (medians of 3) but not
  averaged over instances. Costs depend on m, which varies by instance.
- Proxy sweeps run on CPU (BLAS, 8 threads); statevector references on the
  A100. On a GPU-less machine the statevector column inflates by orders of
  magnitude while the proxy column is unchanged — the comparison is
  machine-dependent in the proxy's favor on commodity hardware.
- GPU single-evaluation timings at small n are latency-bound (~0.2 ms floor),
  so per-evaluation scaling in n is not visible until n ≥ 18.
- sampled_N here is the CPU implementation; it has not been GPU-ported.
- The first Slurm attempt (job 9807239, `log_e013_9807239.log`) crashed at
  n=20 in exact_N: the homodist kernel launched one thread per (x,y) pair and
  4^20/256 blocks overflows CUDA's UInt32 grid limit. The kernel was rewritten
  as a grid-stride loop (commit `ce337950`; `check_gpu_homodist.jl` validates
  GPU vs CPU exactly at n=10, 12, asserted in-job before timing). Job
  9807549 (5:53 elapsed) is the complete run.

## Reproduce

```bash
sbatch research/experiments/013_timing-benchmark/run.sb       # on HPCC (A100)
E13_SMOKE=1 julia --project research/experiments/013_timing-benchmark/run.jl  # local smoke
```

Output: `results.csv` (long format: n, m, op, seconds, reps). Hardware is
printed at the top of the job log (`log_e013_<jobid>.log`).
