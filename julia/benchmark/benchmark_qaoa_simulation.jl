#=
Benchmark QAOA statevector simulation across backends, both Float32 and Float64:
  1. CPU            — qaoa_expectation (pure Julia FUR)
  2. KA-CPU         — gpu_qaoa_expectation with plain Arrays (KernelAbstractions CPU)
  3. GPU            — gpu_qaoa_expectation with device arrays (per-qubit kernel)
  4. GPU-batched    — gpu_qaoa_expectation_batched (shared-memory batched kernel)

All backends are extrapolated using O(2^n · n) scaling once they exceed 250 ms.

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
extrapolate_threshold = 0.25  # seconds; backends are extrapolated once they exceed this
p = 4
edge_prob = 0.5
seed = 42

# ─── Print header ────────────────────────────────────────────────────────────

gpu_label = gpu_name == "none" ? "GPU (—)" : "GPU ($gpu_name)"
batched_label = gpu_name == "none" ? "Batched (—)" : "Batched ($gpu_name)"

TABLE_WIDTH = 180

println("=" ^ TABLE_WIDTH)
println("QAOA Statevector Simulation Benchmark")
println("=" ^ TABLE_WIDTH)
println("  p = $p layers, edge_prob = $edge_prob, seed = $seed")
println("  Backends extrapolated (†) via O(2^n·n) scaling once they exceed $(Int(extrapolate_threshold*1000)) ms")
println("  KA-CPU threads: ", Threads.nthreads())
println("  GPU backend: $gpu_name", gpu_has_f64 ? "" : gpu_name == "none" ? "" : " (Float32 only)")
println("-" ^ TABLE_WIDTH)
@printf("  %4s  %6s  │  %14s  %14s  │  %14s  %14s  │  %14s  %14s  │  %14s  %14s\n",
        "n", "edges",
        "CPU f64", "CPU f32",
        "KA-CPU f64", "KA-CPU f32",
        gpu_has_f64 ? "$gpu_label f64" : "$gpu_label f64 —",
        "$gpu_label f32",
        gpu_has_f64 ? "$batched_label f64" : "$batched_label f64 —",
        "$batched_label f32")
println("-" ^ TABLE_WIDTH)

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
        gpu_qaoa_expectation_batched(gpu_warmup_f32, 4, warmup_γ, warmup_β; group_size=2)
        if gpu_has_f64
            gpu_warmup_f64 = gpu_array_type(warmup_costs_f64)
            gpu_qaoa_expectation(gpu_warmup_f64, 4, warmup_γ, warmup_β)
            gpu_qaoa_expectation_batched(gpu_warmup_f64, 4, warmup_γ, warmup_β; group_size=2)
        end
    end
end

# ─── Benchmark loop ──────────────────────────────────────────────────────────

# Track last measured (n, time) per backend for extrapolation
last_measured = Dict{String, Tuple{Int, Float64}}()

"""
Measure or extrapolate a backend.  Returns (time_seconds, is_extrapolated).
If the backend already exceeded the threshold at a previous n, extrapolate from then.
"""
function measure_or_extrapolate!(key, last_measured, n, f; bm_seconds, bm_samples, threshold=extrapolate_threshold)
    if haskey(last_measured, key)
        n0, t0 = last_measured[key]
        if t0 >= threshold
            return (extrapolate_time(t0, n0, n), true)
        end
    end
    t = f(bm_seconds, bm_samples)
    last_measured[key] = (n, t)
    return (t, false)
end

for n in n_values
    edges = random_edges(n, seed, edge_prob)
    num_edges = length(edges)

    # Random QAOA parameters
    param_rng = MersenneTwister(seed + 1)
    γs = rand(param_rng, p) .* 2π
    βs = rand(param_rng, p) .* 2π

    # Fewer samples for large problems where each run is slow
    bm_samples = n <= 16 ? 100 : 3
    bm_seconds = n <= 16 ? 5   : 1

    costs_f64 = maxcut_costs(n, edges)
    costs_f32 = Float32.(costs_f64)

    # ── CPU benchmarks ───────────────────────────────────────────────────
    t_cpu_f64, ext = measure_or_extrapolate!("cpu_f64", last_measured, n,
        (s, samp) -> @belapsed(qaoa_expectation($costs_f64, $n, $γs, $βs), seconds=s, samples=samp);
        bm_seconds, bm_samples)
    cpu_f64_str = fmt_time(t_cpu_f64; extrapolated=ext)

    t_cpu_f32, ext = measure_or_extrapolate!("cpu_f32", last_measured, n,
        (s, samp) -> @belapsed(qaoa_expectation($costs_f32, $n, $γs, $βs), seconds=s, samples=samp);
        bm_seconds, bm_samples)
    cpu_f32_str = fmt_time(t_cpu_f32; extrapolated=ext)

    # ── KA-CPU benchmarks ────────────────────────────────────────────────
    t_ka_f64, ext = measure_or_extrapolate!("ka_f64", last_measured, n,
        (s, samp) -> @belapsed(gpu_qaoa_expectation($costs_f64, $n, $γs, $βs), seconds=s, samples=samp);
        bm_seconds, bm_samples)
    ka_f64_str = fmt_time(t_ka_f64; extrapolated=ext)

    t_ka_f32, ext = measure_or_extrapolate!("ka_f32", last_measured, n,
        (s, samp) -> @belapsed(gpu_qaoa_expectation($costs_f32, $n, $γs, $βs), seconds=s, samples=samp);
        bm_seconds, bm_samples)
    ka_f32_str = fmt_time(t_ka_f32; extrapolated=ext)

    # ── GPU benchmarks (if available) ────────────────────────────────────
    gpu_f64_str = "—"
    gpu_f32_str = "—"
    bat_f64_str = "—"
    bat_f32_str = "—"

    if gpu_backend !== nothing
        gpu_costs_f32 = gpu_array_type(Float32.(costs_f64))
        t_gpu_f32, ext = measure_or_extrapolate!("gpu_f32", last_measured, n,
            (s, samp) -> @belapsed(gpu_qaoa_expectation($gpu_costs_f32, $n, $γs, $βs), seconds=s, samples=samp);
            bm_seconds, bm_samples)
        gpu_f32_str = fmt_time(t_gpu_f32; extrapolated=ext)

        t_bat_f32, ext = measure_or_extrapolate!("bat_f32", last_measured, n,
            (s, samp) -> @belapsed(gpu_qaoa_expectation_batched($gpu_costs_f32, $n, $γs, $βs), seconds=s, samples=samp);
            bm_seconds, bm_samples)
        bat_f32_str = fmt_time(t_bat_f32; extrapolated=ext)

        if gpu_has_f64
            gpu_costs_f64 = gpu_array_type(costs_f64)
            t_gpu_f64, ext = measure_or_extrapolate!("gpu_f64", last_measured, n,
                (s, samp) -> @belapsed(gpu_qaoa_expectation($gpu_costs_f64, $n, $γs, $βs), seconds=s, samples=samp);
                bm_seconds, bm_samples)
            gpu_f64_str = fmt_time(t_gpu_f64; extrapolated=ext)

            t_bat_f64, ext = measure_or_extrapolate!("bat_f64", last_measured, n,
                (s, samp) -> @belapsed(gpu_qaoa_expectation_batched($gpu_costs_f64, $n, $γs, $βs), seconds=s, samples=samp);
                bm_seconds, bm_samples)
            bat_f64_str = fmt_time(t_bat_f64; extrapolated=ext)
        end
    end

    @printf("  %4d  %6d  │  %14s  %14s  │  %14s  %14s  │  %14s  %14s  │  %14s  %14s\n",
            n, num_edges,
            cpu_f64_str, cpu_f32_str,
            ka_f64_str, ka_f32_str,
            gpu_f64_str, gpu_f32_str,
            bat_f64_str, bat_f32_str)
end

println("=" ^ TABLE_WIDTH)
println("  † = extrapolated from last measured n using O(2^n·n) scaling")
println()
