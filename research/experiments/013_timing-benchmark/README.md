# 013 — Timing benchmark: what does each stage of parameter setting cost?

**Question:** The paper's cost claims are asymptotic ($O(nm^2)$ proxy layers
vs.\ $O(2^n)$ statevector layers vs.\ $O(4^n)$ exact distributions, $O(S m 2^n)$
sampled). What are the actual wall-clock numbers on the hardware the
experiments ran on, so the paper can print one honest timing table?

**Answer: The exact N is the pipeline's bottleneck and the GPU statevector
is launch-bound at these sizes.** exact N costs 21.8 s at n=20 (clean O(4ⁿ):
×13–17 per Δn=2) against the 0.48 s brute-force p=1 ceiling it competes
with; sampled N (S=10) removes that wall at 50× lower cost (0.44 s). A
single GPU statevector pass is ~0.2–0.3 ms *independent of n* up to 20
(kernel-launch-bound), so at simulable sizes exhaustive grid search matches
the proxy pipeline end to end — even at p=20 ramps (est. 16 s brute force
vs 13.4 s proxy sweep + 0.4 s sampled N). The pipeline's cost advantage
opens only as 2ⁿ outruns launch overhead, around n ≈ 22–24 at p=20
(extrapolation). At these sizes the compression's practical value is the
error calculus plus the reuse of one N across every depth and schedule
grid.

Run: job 9807549 (5:53 elapsed). Job 9807239 failed at n=20: the GPU
homodist kernel launched one thread per (x,y) pair and overflowed CUDA's
UInt32 grid limit; fixed by grid-striding the kernel (QOKit commit
`ce337950`), with GPU-vs-CPU equality asserted in-job before timing
(`check_gpu_homodist.jl`, max abs error 0.0 at n=10,12).

## Method

One G(n, 0.5) instance per n ∈ {12, 14, 16, 18, 20} (seed 20260612 + n).
Median wall-clock over 3 runs after a warmup call (compilation excluded;
300 s budget per op). Operations:

| op | what it is |
|----|------------|
| `exact_N` | empirical N(v';d,v), GPU direct homodist, O(4ⁿ) |
| `sampled_N` | stratified sampled N, S=10 per class, CPU, O(S·m·2ⁿ) |
| `proxy_sweep_p1` | full parameter-setting decision over the 40×40 (γ,β) grid: `QAOA_proxy_multi` + `expectation` + argmax |
| `proxy_sweep_p3` | same over the 8⁴ ramp-endpoint grid at p=3 |
| `proxy_sweep_p20` | same 8⁴ schedules at p=20 |
| `sv_single_p1` | one statevector expectation (GPU, costs pre-uploaded) |
| `sv_single_p20` | one 20-layer ramp statevector expectation (GPU) |
| `ceiling_p1` | brute-force grid ceiling = 1600 statevector evaluations |
| `ceiling_p3` | brute-force ramp ceiling = 4096 statevector evaluations |

Grids and bounds identical to experiments 002/004/010 (γ ∈ [0,π], β ∈ [0,π/2]
at p=1; ramp endpoints γ ∈ [0.05,1.6], β ∈ [0.05,0.8]). GPU timings exclude
the one-time cost upload (steady-state regime, matching how the ceilings were
computed in 002/004/010). Proxy sweeps run on CPU (BLAS) — the proxy never
needs the GPU.

## Reproduce

```bash
sbatch research/experiments/013_timing-benchmark/run.sb       # on HPCC
E13_SMOKE=1 julia --project research/experiments/013_timing-benchmark/run.jl  # local smoke
```

Hardware is printed at the top of the job log (`log_e013_<jobid>.log`).
