#=
HPC pilot: validate every GPU code path Phase 3 depends on, before any large
job is built on them. Fail-fast: any disagreement with the CPU reference
aborts with a nonzero exit.

Checks:
  1. CUDA is functional with Float64.
  2. gpu_maxcut_costs == maxcut_costs (n=16).
  3. gpu_qaoa_statevector and gpu_qaoa_statevector_batched == qaoa_statevector
     to 1e-12 (n=16, p=3).
  4. gpu_qaoa_expectation_batched == qaoa_expectation over a small (γ,β) grid.
  5. gpu_get_homogeneous_distribution_from_costs_direct == CPU direct (n=12).
  6. Timing report: CPU vs GPU statevector at n=20 (the Phase-3 working size).

Run via Slurm (scripts/hpc/pilot_gpu.sb) or on a GPU dev node:
  julia scripts/hpc/pilot_gpu.jl        # NOTE: no --project; CUDA loads from
                                        # the default env first (common.jl pattern)
=#

using CUDA, KernelAbstractions
@assert CUDA.functional() "CUDA not functional on this node"
println("GPU: ", CUDA.name(CUDA.device()))

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using JuliaQAOA
using Random: MersenneTwister

elapsed(f) = (t0 = time_ns(); r = f(); ((time_ns() - t0) / 1e9, r))

# Float64 roundtrip (rules out unsupported-precision silent failure)
let v = CUDA.CuArray(Float64[1.0, 2.0, 3.0])
    @assert sum(Array(v)) == 6.0
end
println("[1/6] CUDA functional with Float64 ✓")

rng = MersenneTwister(20260612)
n = 16
edges = erdos_renyi_edges(n, 0.5; rng)
m = length(edges)
costs = maxcut_costs(n, edges)

costs_gpu_computed = Array(gpu_maxcut_costs(n, edges; backend=CUDABackend()))
@assert costs_gpu_computed == costs "gpu_maxcut_costs disagrees with CPU"
println("[2/6] gpu_maxcut_costs agrees (n=$n, m=$m) ✓")

γs = [0.2, 0.5, 0.9]
βs = [0.4, 0.25, 0.1]
ψ_cpu = qaoa_statevector(costs, n, γs, βs)
costs_dev = CUDA.CuArray(costs)
for (label, f) in (("plain", gpu_qaoa_statevector), ("batched", gpu_qaoa_statevector_batched))
    ψ_gpu = Array(f(costs_dev, n, γs, βs))
    err = maximum(abs.(ψ_gpu .- ψ_cpu))
    @assert err < 1e-12 "gpu statevector ($label) max abs error $err > 1e-12"
    println("[3/6] gpu_qaoa_statevector ($label) agrees to ", err, " ✓")
end

worst = 0.0
for γ in (0.1, 0.7, 2.0), β in (0.05, 0.6, 1.3)
    e_cpu = qaoa_expectation(costs, n, [γ], [β])
    e_gpu = gpu_qaoa_expectation_batched(costs_dev, n, [γ], [β])
    global worst = max(worst, abs(e_cpu - e_gpu))
end
@assert worst < 1e-9 "gpu expectation max abs error $worst > 1e-9"
println("[4/6] gpu_qaoa_expectation_batched agrees to ", worst, " ✓")

n12 = 12
edges12 = erdos_renyi_edges(n12, 0.5; rng)
m12 = length(edges12)
costs12 = maxcut_costs(n12, edges12)
N_cpu = get_homogeneous_distribution_from_costs_direct(costs12, m12, n12)
N_gpu = Array(gpu_get_homogeneous_distribution_from_costs_direct(costs12, m12, n12))
err_N = maximum(abs.(N_cpu .- N_gpu))
@assert err_N < 1e-10 "gpu homodist max abs error $err_N"
println("[5/6] gpu homodist (direct) agrees to ", err_N, " ✓")

# Timing at the Phase-3 working size
n20 = 20
edges20 = erdos_renyi_edges(n20, 0.5; rng)
costs20 = maxcut_costs(n20, edges20)
costs20_dev = CUDA.CuArray(costs20)
γs20, βs20 = linear_ramp(0.1, 0.8, 0.6, 0.1, 10)
gpu_qaoa_expectation_batched(costs20_dev, n20, γs20, βs20)   # warm up compile
t_cpu, e_cpu = elapsed(() -> qaoa_expectation(costs20, n20, γs20, βs20))
t_gpu, e_gpu = elapsed(() -> gpu_qaoa_expectation_batched(costs20_dev, n20, γs20, βs20))
@assert abs(e_cpu - e_gpu) < 1e-8
println("[6/6] n=$n20 p=10 expectation: CPU ", round(t_cpu, digits=3), "s,  GPU ",
        round(t_gpu, digits=4), "s  (", round(t_cpu / t_gpu, digits=1), "× speedup) ✓")

println("\nPILOT PASSED — all GPU paths agree with CPU references.")
