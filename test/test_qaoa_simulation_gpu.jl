using JuliaQAOA, Test
using KernelAbstractions
import Random: MersenneTwister

# ─── GPU auto-detection ────────────────────────────────────────────────────
# `using` inside try blocks needs `@eval` to bring the module into scope,
# and `global` to assign to file-level variables from soft scope.

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

# Use Float32 for backends without native Float64 support
gpu_float_type = (gpu_name in ("oneAPI", "Metal")) ? Float32 : Float64

@info "GPU backend detected: $gpu_name"

# ─── Helper: generate random ER graph edges ─────────────────────────────────
function random_edges(n, seed)
    rng = MersenneTwister(seed)
    edges = Tuple{Int,Int}[]
    for i in 0:(n-2), j in (i+1):(n-1)
        if rand(rng) < 0.5
            push!(edges, (i, j))
        end
    end
    return edges
end

# ─── Tests that require a GPU ──────────────────────────────────────────────

@testset "GPU QAOA simulation matches CPU ($gpu_name Float64)" begin
    if gpu_backend === nothing || gpu_float_type != Float64
        @warn "Skipping Float64 GPU test (backend: $gpu_name)"
        return
    end

    for seed in 0:4
        edges = random_edges(10, seed)
        p = 3
        param_rng = MersenneTwister(seed + 100)
        γs = rand(param_rng, p) .* 2π
        βs = rand(param_rng, p) .* 2π

        cpu_costs = maxcut_costs(10, edges)
        cpu_state = qaoa_statevector(cpu_costs, 10, γs, βs)
        cpu_exp = qaoa_expectation(cpu_costs, 10, γs, βs)

        gpu_costs = gpu_array_type(cpu_costs)
        gpu_state = gpu_qaoa_statevector(gpu_costs, 10, γs, βs)
        gpu_exp = gpu_qaoa_expectation(gpu_costs, 10, γs, βs)

        @test Array(gpu_state) ≈ cpu_state atol=1e-10
        @test gpu_exp ≈ cpu_exp atol=1e-10
    end
end

@testset "GPU QAOA simulation ($gpu_name Float32)" begin
    if gpu_backend === nothing
        @warn "No GPU available, skipping"
        return
    end

    for seed in 0:4
        edges = random_edges(10, seed)
        p = 3
        param_rng = MersenneTwister(seed + 100)
        γs = rand(param_rng, p) .* 2π
        βs = rand(param_rng, p) .* 2π

        costs_f32 = Float32.(maxcut_costs(10, edges))
        cpu_state = qaoa_statevector(costs_f32, 10, γs, βs)
        cpu_exp = qaoa_expectation(costs_f32, 10, γs, βs)

        gpu_costs = gpu_array_type(costs_f32)
        gpu_state = gpu_qaoa_statevector(gpu_costs, 10, γs, βs)
        gpu_exp = gpu_qaoa_expectation(gpu_costs, 10, γs, βs)

        @test Array(gpu_state) ≈ cpu_state atol=1e-4
        @test gpu_exp ≈ cpu_exp atol=1e-4
    end
end

@testset "GPU phase gate matches CPU ($gpu_name)" begin
    if gpu_backend === nothing
        @warn "No GPU available, skipping"
        return
    end

    edges = random_edges(8, 42)
    T = gpu_float_type
    CT = Complex{T}
    costs = T.(maxcut_costs(8, edges))
    γ = 1.23

    cpu_state = fill(CT(1 / sqrt(T(1 << 8))), 1 << 8)
    apply_phase_gate!(cpu_state, costs, γ)

    gpu_state = gpu_array_type(fill(CT(1 / sqrt(T(1 << 8))), 1 << 8))
    gpu_apply_phase_gate!(gpu_state, gpu_array_type(costs), γ)

    atol = T == Float32 ? 1e-5 : 1e-12
    @test Array(gpu_state) ≈ cpu_state atol=atol
end

@testset "GPU X-mixer matches CPU ($gpu_name)" begin
    if gpu_backend === nothing
        @warn "No GPU available, skipping"
        return
    end

    n = 8
    β = 0.789
    T = gpu_float_type
    CT = Complex{T}

    rng = MersenneTwister(99)
    cpu_state = CT.(randn(rng, ComplexF64, 1 << n))
    cpu_state ./= sqrt(T(sum(abs2, cpu_state)))

    gpu_state = gpu_array_type(copy(cpu_state))

    apply_x_mixer!(cpu_state, β, n)
    gpu_apply_x_mixer!(gpu_state, β, n)

    atol = T == Float32 ? 1e-4 : 1e-12
    @test Array(gpu_state) ≈ cpu_state atol=atol
end

@testset "gpu_maxcut_costs convenience function ($gpu_name)" begin
    if gpu_backend === nothing
        @warn "No GPU available, skipping"
        return
    end

    edges = [(0,1), (1,2), (2,3), (0,3), (1,3)]
    cpu_costs = maxcut_costs(6, edges)

    if gpu_float_type == Float64
        gpu_costs64 = gpu_maxcut_costs(6, edges; backend=gpu_backend)
        @test Array(gpu_costs64) == cpu_costs
        @test eltype(gpu_costs64) == Float64
    end

    gpu_costs32 = gpu_maxcut_costs(6, edges; backend=gpu_backend, T=Float32)
    @test Array(gpu_costs32) ≈ Float32.(cpu_costs)
    @test eltype(gpu_costs32) == Float32
end

# ─── Tests using CPU backend (always run) ────────────────────────────────────
# KernelAbstractions runs on plain Arrays via CPU backend, so we can validate
# kernel correctness without a GPU.

@testset "gpu_maxcut_costs on CPU backend" begin
    using KernelAbstractions: CPU

    for seed in 0:4
        edges = random_edges(10, seed)
        cpu_costs = maxcut_costs(10, edges)

        ka_costs64 = gpu_maxcut_costs(10, edges; backend=CPU(), T=Float64)
        @test ka_costs64 == cpu_costs
        @test eltype(ka_costs64) == Float64

        ka_costs32 = gpu_maxcut_costs(10, edges; backend=CPU(), T=Float32)
        @test ka_costs32 == Float32.(cpu_costs)
        @test eltype(ka_costs32) == Float32
    end
end

@testset "KA kernels on CPU backend (Float64)" begin
    for seed in 0:4
        edges = random_edges(10, seed)
        p = 3
        param_rng = MersenneTwister(seed + 100)
        γs = rand(param_rng, p) .* 2π
        βs = rand(param_rng, p) .* 2π

        costs = maxcut_costs(10, edges)
        cpu_state = qaoa_statevector(costs, 10, γs, βs)
        cpu_exp = qaoa_expectation(costs, 10, γs, βs)

        ka_state = gpu_qaoa_statevector(costs, 10, γs, βs)
        ka_exp = gpu_qaoa_expectation(costs, 10, γs, βs)

        @test ka_state ≈ cpu_state atol=1e-12
        @test ka_exp ≈ cpu_exp atol=1e-10
    end
end

@testset "KA batched X-mixer on CPU backend (Float64)" begin
    for seed in 0:4
        edges = random_edges(10, seed)
        p = 3
        param_rng = MersenneTwister(seed + 100)
        γs = rand(param_rng, p) .* 2π
        βs = rand(param_rng, p) .* 2π

        costs = maxcut_costs(10, edges)
        cpu_exp = qaoa_expectation(costs, 10, γs, βs)

        for gs in [2, 5, 10]
            bat_exp = gpu_qaoa_expectation_batched(costs, 10, γs, βs; group_size=gs)
            @test bat_exp ≈ cpu_exp atol=1e-10
        end
    end
end

@testset "KA batched X-mixer on CPU backend (Float32)" begin
    for seed in 0:4
        edges = random_edges(10, seed)
        p = 3
        param_rng = MersenneTwister(seed + 100)
        γs = rand(param_rng, p) .* 2π
        βs = rand(param_rng, p) .* 2π

        costs32 = Float32.(maxcut_costs(10, edges))
        cpu_exp = qaoa_expectation(costs32, 10, γs, βs)

        for gs in [2, 5, 10]
            bat_exp = gpu_qaoa_expectation_batched(costs32, 10, γs, βs; group_size=gs)
            @test bat_exp ≈ cpu_exp atol=1e-4
        end
    end
end

@testset "GPU batched QAOA matches CPU ($gpu_name Float64)" begin
    if gpu_backend === nothing || gpu_float_type != Float64
        @warn "Skipping Float64 GPU batched test (backend: $gpu_name)"
        return
    end

    for seed in 0:4
        edges = random_edges(10, seed)
        p = 3
        param_rng = MersenneTwister(seed + 100)
        γs = rand(param_rng, p) .* 2π
        βs = rand(param_rng, p) .* 2π

        cpu_costs = maxcut_costs(10, edges)
        cpu_exp = qaoa_expectation(cpu_costs, 10, γs, βs)

        gpu_costs = gpu_array_type(cpu_costs)
        for gs in [2, 5, 10]
            gpu_exp = gpu_qaoa_expectation_batched(gpu_costs, 10, γs, βs; group_size=gs)
            @test gpu_exp ≈ cpu_exp atol=1e-10
        end
    end
end

@testset "GPU batched QAOA ($gpu_name Float32)" begin
    if gpu_backend === nothing
        @warn "No GPU available, skipping"
        return
    end

    for seed in 0:4
        edges = random_edges(10, seed)
        p = 3
        param_rng = MersenneTwister(seed + 100)
        γs = rand(param_rng, p) .* 2π
        βs = rand(param_rng, p) .* 2π

        costs_f32 = Float32.(maxcut_costs(10, edges))
        cpu_exp = qaoa_expectation(costs_f32, 10, γs, βs)

        gpu_costs = gpu_array_type(costs_f32)
        for gs in [2, 5, 10]
            gpu_exp = gpu_qaoa_expectation_batched(gpu_costs, 10, γs, βs; group_size=gs)
            @test gpu_exp ≈ cpu_exp atol=1e-4
        end
    end
end

@testset "KA kernels on CPU backend (Float32)" begin
    for seed in 0:4
        edges = random_edges(10, seed)
        p = 3
        param_rng = MersenneTwister(seed + 100)
        γs = rand(param_rng, p) .* 2π
        βs = rand(param_rng, p) .* 2π

        costs32 = Float32.(maxcut_costs(10, edges))
        cpu_state = qaoa_statevector(costs32, 10, γs, βs)
        cpu_exp = qaoa_expectation(costs32, 10, γs, βs)

        ka_state = gpu_qaoa_statevector(costs32, 10, γs, βs)
        ka_exp = gpu_qaoa_expectation(costs32, 10, γs, βs)

        @test eltype(ka_state) == ComplexF32
        @test ka_state ≈ cpu_state atol=1e-5
        @test ka_exp ≈ cpu_exp atol=1e-4
    end
end
