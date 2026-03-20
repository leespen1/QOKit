#=
Benchmark Julia vs Python QAOA statevector simulation.

Compares:
  1. Julia CPU        — qaoa_expectation (pure Julia FUR)
  2. Julia KA-CPU     — gpu_qaoa_expectation with Arrays (KernelAbstractions CPU)
  3. Julia GPU        — gpu_qaoa_expectation on device (auto-detected)
  4. Python (C)       — QOKit C-compiled FUR backend
  5. Python (GPU)     — QOKit CUDA/Numba FUR backend
  6. Python (pure)    — QOKit pure-Python FUR backend

Usage:
    julia --project=. benchmark/benchmark_julia_vs_python.jl
=#

include(joinpath(@__DIR__, "..", "find_python.jl"))

using JuliaQAOA
using BenchmarkTools
using Printf
using PythonCall
import Random: MersenneTwister

# ── Python setup ──────────────────────────────────────────────────────────
const np = pyimport("numpy")
const nx = pyimport("networkx")
const mc = pyimport("qokit.maxcut")
const qaoa_sim = pyimport("grips.QAOA_simulator")
const py_fur = pyimport("qokit.fur")
const py_time = pyimport("time")

py_available = pyconvert(Vector{String}, py_fur.get_available_simulator_names())
has_py_c = "c" in py_available
has_py_gpu = "gpu" in py_available

function make_nx_graph(n::Int, edges::Vector{Tuple{Int,Int}})
    G = nx.Graph()
    G.add_nodes_from(pylist(0:(n-1)))
    for (i, j) in edges
        G.add_edge(i, j)
    end
    return G
end

"""
Time a Python QAOA expectation call by running it `repeats` times and
returning the median wall-clock time (in seconds).
"""
function time_python_expectation(n, terms, γs, βs, simulator_name; repeats=5)
    py_gamma = np.array(γs)
    py_beta = np.array(βs)

    # Create simulator once (amortize setup cost)
    sim = qaoa_sim.get_simulator(n, terms, simulator_name=simulator_name)

    # Warmup
    qaoa_sim.get_expectation(n, terms, py_gamma, py_beta,
                             sim=sim, simulator_name=simulator_name)

    times = Float64[]
    for _ in 1:repeats
        t0 = pyconvert(Float64, py_time.perf_counter())
        qaoa_sim.get_expectation(n, terms, py_gamma, py_beta,
                                 sim=sim, simulator_name=simulator_name)
        t1 = pyconvert(Float64, py_time.perf_counter())
        push!(times, t1 - t0)
    end
    return sort(times)[div(repeats, 2) + 1]  # median
end

# ─── GPU auto-detection ────────────────────────────────────────────────────

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
gpu_T = gpu_has_f64 ? Float64 : Float32

# ─── Configuration ─────────────────────────────────────────────────────────

n_values = [8, 12, 16, 20, 24, 28]
extrapolate_threshold = 0.25  # seconds; backends are extrapolated once they exceed this
p = 4
edge_prob = 0.5
seed = 42

# ─── Helpers ──────────────────────────────────────────────────────────────

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
"""
function extrapolate_time(t_measured, n_measured, n_target)
    return t_measured * (exp2(n_target) * n_target) / (exp2(n_measured) * n_measured)
end

# ─── Print header ──────────────────────────────────────────────────────────

# All columns use Float64 except GPU when the backend lacks support
gpu_col = if gpu_name == "none"
    "Jl GPU (—)"
elseif !gpu_has_f64
    "Jl GPU f32"
else
    "Jl GPU"
end

println("=" ^ 115)
println("QAOA Statevector Simulation: Julia vs Python (all Float64)")
println("=" ^ 115)
println("  p = $p layers, edge_prob = $edge_prob, seed = $seed")
println("  Backends extrapolated (†) via O(2^n·n) scaling once they exceed $(Int(extrapolate_threshold*1000)) ms")
println("  Julia KA-CPU threads: ", Threads.nthreads())
println("  Julia GPU backend: $gpu_name", !gpu_has_f64 && gpu_name != "none" ? " (Float32 — no native Float64)" : "")
println("  Python backends: ", join(py_available, ", "), " (always Float64)")
println("-" ^ 115)

@printf("  %4s  %6s  │  %14s  %14s  %14s  │  %14s  %14s  %14s\n",
        "n", "edges", "Jl CPU", "Jl KA-CPU", gpu_col,
        has_py_c ? "Py C" : "Py C (—)",
        has_py_gpu ? "Py GPU" : "Py GPU (—)",
        "Py pure")
println("-" ^ 115)

# ─── JIT warmup (small n to trigger compilation without large allocations) ──

let
    warmup_costs = maxcut_costs(4, [(0,1),(1,2),(2,3)])
    warmup_γ = [1.0]
    warmup_β = [1.0]
    qaoa_expectation(warmup_costs, 4, warmup_γ, warmup_β)
    gpu_qaoa_expectation(warmup_costs, 4, warmup_γ, warmup_β)
    if gpu_backend !== nothing
        gpu_warmup = gpu_array_type(gpu_T.(warmup_costs))
        gpu_qaoa_expectation(gpu_warmup, 4, warmup_γ, warmup_β)
    end
end

# ─── Benchmark loop ────────────────────────────────────────────────────────

# Track last measured (n, time) per non-GPU backend for extrapolation
last_measured = Dict{String, Tuple{Int, Float64}}()

"""
Measure or extrapolate a non-GPU backend.  Returns (time_seconds, is_extrapolated).
If the backend already exceeded the threshold at a previous n, extrapolate from then.
"""
function measure_or_extrapolate!(key, last_measured, n, f; bm_seconds, bm_samples, threshold=extrapolate_threshold)
    if haskey(last_measured, key)
        n0, t0 = last_measured[key]
        if t0 >= threshold
            # Already exceeded threshold — extrapolate
            return (extrapolate_time(t0, n0, n), true)
        end
    end
    t = f(bm_seconds, bm_samples)
    last_measured[key] = (n, t)
    return (t, false)
end

function measure_or_extrapolate_py!(key, last_measured, n, terms, γs, βs, sim_name; py_repeats, threshold=extrapolate_threshold)
    if haskey(last_measured, key)
        n0, t0 = last_measured[key]
        if t0 >= threshold
            return (extrapolate_time(t0, n0, n), true)
        end
    end
    t = time_python_expectation(n, terms, γs, βs, sim_name; repeats=py_repeats)
    last_measured[key] = (n, t)
    return (t, false)
end

for n in n_values
    # Generate random graph
    rng = MersenneTwister(seed)
    edges = Tuple{Int,Int}[]
    for i in 0:(n-2), j in (i+1):(n-1)
        if rand(rng) < edge_prob
            push!(edges, (i, j))
        end
    end
    num_edges = length(edges)

    # Random QAOA parameters
    param_rng = MersenneTwister(seed + 1)
    γs = rand(param_rng, p) .* 2π
    βs = rand(param_rng, p) .* 2π

    # Python graph and terms (created once per n)
    G = make_nx_graph(n, edges)
    terms = mc.get_maxcut_terms(G)

    # Fewer samples for large problems where each run is slow
    bm_samples = n <= 16 ? 100 : 3
    bm_seconds = n <= 16 ? 5   : 1
    py_repeats = n <= 16 ? 5   : 3

    cpu_costs = maxcut_costs(n, edges)

    # ── Julia CPU (Float64) ───────────────────────────────────────────
    t_jl_cpu, ext_cpu = measure_or_extrapolate!("jl_cpu", last_measured, n,
        (s, samp) -> @belapsed(qaoa_expectation($cpu_costs, $n, $γs, $βs), seconds=s, samples=samp);
        bm_seconds, bm_samples)
    jl_cpu_str = fmt_time(t_jl_cpu; extrapolated=ext_cpu)

    # ── Julia KA-CPU (Float64) ────────────────────────────────────────
    t_jl_ka, ext_ka = measure_or_extrapolate!("jl_ka", last_measured, n,
        (s, samp) -> @belapsed(gpu_qaoa_expectation($cpu_costs, $n, $γs, $βs), seconds=s, samples=samp);
        bm_seconds, bm_samples)
    jl_ka_str = fmt_time(t_jl_ka; extrapolated=ext_ka)

    # ── Julia GPU (Float64 if supported, else Float32) ────────────────
    jl_gpu_str = "—"
    if gpu_backend !== nothing
        gpu_costs = gpu_array_type(gpu_T.(cpu_costs))
        t_jl_gpu, ext_gpu = measure_or_extrapolate!("jl_gpu", last_measured, n,
            (s, samp) -> @belapsed(gpu_qaoa_expectation($gpu_costs, $n, $γs, $βs), seconds=s, samples=samp);
            bm_seconds, bm_samples)
        jl_gpu_str = fmt_time(t_jl_gpu; extrapolated=ext_gpu)
    end

    # ── Python C backend (Float64) ────────────────────────────────────
    py_c_str = "—"
    if has_py_c
        t_py_c, ext_py_c = measure_or_extrapolate_py!("py_c", last_measured, n, terms, γs, βs, "c";
            py_repeats)
        py_c_str = fmt_time(t_py_c; extrapolated=ext_py_c)
    end

    # ── Python GPU backend (Float64) ──────────────────────────────────
    py_gpu_str = "—"
    if has_py_gpu
        try
            t_py_gpu, ext_py_gpu = measure_or_extrapolate_py!("py_gpu", last_measured, n, terms, γs, βs, "gpu";
                py_repeats)
            py_gpu_str = fmt_time(t_py_gpu; extrapolated=ext_py_gpu)
        catch e
            @warn "Python GPU backend failed, skipping for remaining sizes" exception=(e, catch_backtrace())
            has_py_gpu = false
        end
    end

    # ── Python pure backend (Float64) ─────────────────────────────────
    t_py_pure, ext_py_pure = measure_or_extrapolate_py!("py_pure", last_measured, n, terms, γs, βs, "python";
        py_repeats)
    py_pure_str = fmt_time(t_py_pure; extrapolated=ext_py_pure)

    @printf("  %4d  %6d  │  %14s  %14s  %14s  │  %14s  %14s  %14s\n",
            n, num_edges,
            jl_cpu_str, jl_ka_str, jl_gpu_str,
            py_c_str, py_gpu_str, py_pure_str)
end

println("=" ^ 115)
println("  † = extrapolated from last measured n using O(2^n·n) scaling")

# ─── Correctness verification ──────────────────────────────────────────────

println("\nCorrectness check (n=12, p=$p):")
rng = MersenneTwister(seed)
n_check = 12
edges_check = Tuple{Int,Int}[]
for i in 0:(n_check-2), j in (i+1):(n_check-1)
    if rand(rng) < edge_prob
        push!(edges_check, (i, j))
    end
end
param_rng = MersenneTwister(seed + 1)
γs_check = rand(param_rng, p) .* 2π
βs_check = rand(param_rng, p) .* 2π

costs_check = maxcut_costs(n_check, edges_check)
jl_exp = qaoa_expectation(costs_check, n_check, γs_check, βs_check)

G_check = make_nx_graph(n_check, edges_check)
terms_check = mc.get_maxcut_terms(G_check)
py_exp = pyconvert(Float64, qaoa_sim.get_expectation(
    n_check, terms_check,
    np.array(γs_check), np.array(βs_check),
    simulator_name="python"
))

@printf("  Julia CPU:    %.12f\n", jl_exp)
@printf("  Python pure:  %.12f\n", py_exp)
@printf("  Difference:   %.2e\n", abs(jl_exp - py_exp))
