#=
Benchmark QAOA statevector simulation across three backends, both Float32 and Float64:
  1. CPU          — qaoa_expectation (pure Julia FUR, Float64 only)
  2. KA-CPU       — gpu_qaoa_expectation with plain Arrays (KernelAbstractions CPU backend)
  3. GPU          — gpu_qaoa_expectation with device arrays (auto-detected backend)

CPU and KA-CPU are benchmarked up to n=20; beyond that, times are extrapolated
using the O(2^n · n) scaling measured from the last two real data points.
GPU is benchmarked for real at all n values up to n=28.

Float64 GPU is skipped on backends that lack native support (oneAPI, Metal).

Usage:
    julia --project=. benchmark/benchmark_qaoa_simulation.jl
=#

using JuliaQAOA
using KernelAbstractions
using BenchmarkTools
using Printf
import Random: MersenneTwister

# ─── GPU auto-detection ──────────────────────────────────────────────────────

gpu_backend = nothing
gpu_array_type = nothing
gpu_name = "none"

try
    @eval using CUDA
    if CUDA.functional()
        global gpu_backend = CUDA.CUDABackend()
        global gpu_array_type = CUDA.CuArray
        global gpu_name = "CUDA"
    end
catch; end

if gpu_backend === nothing
    try
        @eval using AMDGPU
        if AMDGPU.functional()
            global gpu_backend = AMDGPU.ROCBackend()
            global gpu_array_type = AMDGPU.ROCArray
            global gpu_name = "ROCm"
        end
    catch; end
end

if gpu_backend === nothing
    try
        @eval using oneAPI
        if oneAPI.functional()
            global gpu_backend = oneAPI.oneAPIBackend()
            global gpu_array_type = oneAPI.oneArray
            global gpu_name = "oneAPI"
        end
    catch; end
end

if gpu_backend === nothing
    try
        @eval using Metal
        if Metal.functional()
            global gpu_backend = Metal.MetalBackend()
            global gpu_array_type = Metal.MtlArray
            global gpu_name = "Metal"
        end
    catch; end
end

# Whether this GPU backend supports Float64
gpu_has_f64 = gpu_name ∉ ("oneAPI", "Metal", "none")

# ─── Configuration ───────────────────────────────────────────────────────────

n_values = [8, 12, 16, 20, 24, 28]
max_cpu_n = 20       # CPU/KA-CPU are measured up to this n; beyond is extrapolated
p = 4
edge_prob = 0.5
seed = 42

# ─── Print header ────────────────────────────────────────────────────────────

gpu_label = gpu_name == "none" ? "GPU (—)" : "GPU ($gpu_name)"

println("=" ^ 130)
println("QAOA Statevector Simulation Benchmark")
println("=" ^ 130)
println("  p = $p layers, edge_prob = $edge_prob, seed = $seed")
println("  CPU/KA-CPU measured up to n=$max_cpu_n; extrapolated (†) beyond via O(2^n·n) scaling")
println("  KA-CPU threads: ", Threads.nthreads())
println("  GPU backend: $gpu_name", gpu_has_f64 ? "" : gpu_name == "none" ? "" : " (Float32 only)")
println("-" ^ 130)
@printf("  %4s  %6s  │  %14s  %14s  │  %14s  %14s  │  %14s  %14s\n",
        "n", "edges",
        "CPU f64", "CPU f32",
        "KA-CPU f64", "KA-CPU f32",
        gpu_has_f64 ? "$gpu_label f64" : "$gpu_label f64 —",
        "$gpu_label f32")
println("-" ^ 130)

# ─── Helpers ─────────────────────────────────────────────────────────────────

function fmt_time(t; extrapolated=false)
    suffix = extrapolated ? " †" : "  "
    if t < 1.0
        return @sprintf("%10.3f ms%s", t * 1000, suffix)
    else
        return @sprintf("%10.3f s %s", t, suffix)
    end
end

"""
Extrapolate a CPU time from n₀ to n₁ using the known O(2^n · n) scaling.
The constant factor cancels: t(n₁)/t(n₀) = (2^n₁ · n₁) / (2^n₀ · n₀).
"""
function extrapolate_time(t_measured, n_measured, n_target)
    return t_measured * (exp2(n_target) * n_target) / (exp2(n_measured) * n_measured)
end

function random_edges(n, seed, edge_prob)
    rng = MersenneTwister(seed)
    edges = Tuple{Int,Int}[]
    for i in 0:(n-2), j in (i+1):(n-1)
        if rand(rng) < edge_prob
            push!(edges, (i, j))
        end
    end
    return edges
end

# ─── JIT warmup (small n to trigger compilation without large allocations) ──

let
    warmup_costs_f64 = maxcut_costs(4, [(0,1),(1,2),(2,3)])
    warmup_costs_f32 = Float32.(warmup_costs_f64)
    warmup_γ = [1.0]
    warmup_β = [1.0]
    qaoa_expectation(warmup_costs_f64, 4, warmup_γ, warmup_β)
    qaoa_expectation(warmup_costs_f32, 4, warmup_γ, warmup_β)
    gpu_qaoa_expectation(warmup_costs_f64, 4, warmup_γ, warmup_β)
    gpu_qaoa_expectation(warmup_costs_f32, 4, warmup_γ, warmup_β)
    if gpu_backend !== nothing
        gpu_warmup_f32 = gpu_array_type(warmup_costs_f32)
        gpu_qaoa_expectation(gpu_warmup_f32, 4, warmup_γ, warmup_β)
        if gpu_has_f64
            gpu_warmup_f64 = gpu_array_type(warmup_costs_f64)
            gpu_qaoa_expectation(gpu_warmup_f64, 4, warmup_γ, warmup_β)
        end
    end
end

# ─── Benchmark loop ──────────────────────────────────────────────────────────

# Store the last measured CPU/KA-CPU times and n for extrapolation
last_cpu = Dict{String, Tuple{Int, Float64}}()  # key → (n, time)

for n in n_values
    edges = random_edges(n, seed, edge_prob)
    num_edges = length(edges)

    # Random QAOA parameters
    param_rng = MersenneTwister(seed + 1)
    γs = rand(param_rng, p) .* 2π
    βs = rand(param_rng, p) .* 2π

    is_cpu_measured = n <= max_cpu_n

    # Fewer samples for large problems where each run is slow
    bm_samples = n <= 16 ? 100 : 3
    bm_seconds = n <= 16 ? 5   : 1

    # ── CPU and KA-CPU benchmarks ─────────────────────────────────────────
    if is_cpu_measured
        costs_f64 = maxcut_costs(n, edges)
        costs_f32 = Float32.(costs_f64)

        t_cpu_f64 = @belapsed qaoa_expectation($costs_f64, $n, $γs, $βs) seconds=bm_seconds samples=bm_samples
        t_cpu_f32 = @belapsed qaoa_expectation($costs_f32, $n, $γs, $βs) seconds=bm_seconds samples=bm_samples
        t_ka_f64  = @belapsed gpu_qaoa_expectation($costs_f64, $n, $γs, $βs) seconds=bm_seconds samples=bm_samples
        t_ka_f32  = @belapsed gpu_qaoa_expectation($costs_f32, $n, $γs, $βs) seconds=bm_seconds samples=bm_samples

        last_cpu["cpu_f64"] = (n, t_cpu_f64)
        last_cpu["cpu_f32"] = (n, t_cpu_f32)
        last_cpu["ka_f64"]  = (n, t_ka_f64)
        last_cpu["ka_f32"]  = (n, t_ka_f32)

        cpu_f64_str = fmt_time(t_cpu_f64)
        cpu_f32_str = fmt_time(t_cpu_f32)
        ka_f64_str  = fmt_time(t_ka_f64)
        ka_f32_str  = fmt_time(t_ka_f32)
    else
        # Extrapolate from last measured point
        n0, t0 = last_cpu["cpu_f64"];  cpu_f64_str = fmt_time(extrapolate_time(t0, n0, n), extrapolated=true)
        n0, t0 = last_cpu["cpu_f32"];  cpu_f32_str = fmt_time(extrapolate_time(t0, n0, n), extrapolated=true)
        n0, t0 = last_cpu["ka_f64"];   ka_f64_str  = fmt_time(extrapolate_time(t0, n0, n), extrapolated=true)
        n0, t0 = last_cpu["ka_f32"];   ka_f32_str  = fmt_time(extrapolate_time(t0, n0, n), extrapolated=true)
    end

    # ── GPU benchmark (always real, if available) ─────────────────────────
    gpu_f64_str = "—"
    gpu_f32_str = "—"

    if gpu_backend !== nothing
        gpu_costs_f32 = gpu_array_type(Float32.(maxcut_costs(n, edges)))
        t_gpu_f32 = @belapsed gpu_qaoa_expectation($gpu_costs_f32, $n, $γs, $βs) seconds=bm_seconds samples=bm_samples
        gpu_f32_str = fmt_time(t_gpu_f32)

        if gpu_has_f64
            gpu_costs_f64 = gpu_array_type(Float64.(maxcut_costs(n, edges)))
            t_gpu_f64 = @belapsed gpu_qaoa_expectation($gpu_costs_f64, $n, $γs, $βs) seconds=bm_seconds samples=bm_samples
            gpu_f64_str = fmt_time(t_gpu_f64)
        end
    end

    @printf("  %4d  %6d  │  %14s  %14s  │  %14s  %14s  │  %14s  %14s\n",
            n, num_edges,
            cpu_f64_str, cpu_f32_str,
            ka_f64_str, ka_f32_str,
            gpu_f64_str, gpu_f32_str)
end

println("=" ^ 130)
println("  † = extrapolated from n=$max_cpu_n using O(2^n · n) scaling")
println()
