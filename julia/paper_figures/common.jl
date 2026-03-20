#=
common.jl — Shared utilities for paper figure reproduction scripts.

Provides:
  - Erdős-Rényi graph generation (as edge lists)
  - JuliaQAOA module loading (which provides maxcut_costs, qaoa_statevector,
    qaoa_expectation, etc.)
  - Plotting utilities
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "find_python.jl"))

using JuliaQAOA
using CairoMakie
using Random
using LinearAlgebra

#==============================================================================#
#                          GPU Backend Selection                                #
#==============================================================================#
# Try CUDA first, then AMDGPU. Require Float64 support (rules out Intel iGPU
# and Apple Metal). Falls back to CPU if no suitable GPU is found.
# KernelAbstractions is loaded alongside whichever GPU backend activates, which
# triggers JuliaQAOAKernelAbstractionsExt and makes gpu_qaoa_* available.

const USE_GPU, _GPU_BACKEND = let
    result = (false, :none)

    # --- Try CUDA ---
    try
        @eval using CUDA, KernelAbstractions
        if CUDA.functional()
            try
                test = CUDA.CuArray(Float64[1.0, 2.0, 3.0])
                _ = sum(Array(test))   # Array() syncs and validates Float64 round-trip
                @info "GPU backend: CUDA (Float64)"
                result = (true, :cuda)
            catch
                @info "CUDA device found but Float64 not supported; trying AMDGPU..."
            end
        end
    catch
    end

    # --- Try AMDGPU (if CUDA unavailable or Float64-incapable) ---
    if !result[1]
        try
            @eval using AMDGPU, KernelAbstractions
            if AMDGPU.functional()
                try
                    test = AMDGPU.ROCArray(Float64[1.0, 2.0, 3.0])
                    _ = sum(Array(test))
                    @info "GPU backend: AMDGPU (Float64)"
                    result = (true, :amdgpu)
                catch
                    @info "AMDGPU device found but Float64 not supported; using CPU"
                end
            end
        catch
        end
    end

    !result[1] && @info "No GPU with Float64 support found; using CPU"
    result
end

"""Transfer a CPU Float64 vector to the active GPU array type, or return as-is for CPU."""
function _to_gpu(v::Vector{Float64})
    _GPU_BACKEND === :cuda   && return CUDA.CuArray(v)
    _GPU_BACKEND === :amdgpu && return AMDGPU.ROCArray(v)
    return v
end

"""QAOA expectation value on the preferred backend (GPU if available, CPU otherwise)."""
function qaoa_expectation_device(costs::Vector{Float64}, n::Int,
                                  γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real})
    USE_GPU || return qaoa_expectation(costs, n, γs, βs)
    return Float64(gpu_qaoa_expectation(_to_gpu(costs), n, γs, βs))
end

"""QAOA statevectors at all layers [0, 1, …, p] on the preferred backend.
Returns a Vector of CPU ComplexF64 statevectors (layer 0 = initial uniform state).
"""
function qaoa_statevector_intermediates_device(costs::Vector{Float64}, n::Int,
                                               γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real})
    p = length(γs)
    initial = fill(ComplexF64(1.0 / sqrt(1 << n)), 1 << n)
    if USE_GPU
        costs_gpu = _to_gpu(costs)
        state = similar(costs_gpu, Complex{eltype(costs_gpu)})
        fill!(state, eltype(state)(1.0 / sqrt(1 << n)))
        states = Vector{ComplexF64}[copy(initial)]
        for ℓ in 1:p
            gpu_apply_phase_gate!(state, costs_gpu, γs[ℓ])
            gpu_apply_x_mixer!(state, βs[ℓ], n)
            push!(states, Vector{ComplexF64}(Array(state)))
        end
        return states
    else
        intermediates = qaoa_statevector(costs, n, γs, βs; return_intermediates=true)
        pushfirst!(intermediates, initial)
        return intermediates
    end
end

#==============================================================================#
#                          Graph Generation                                     #
#==============================================================================#

"""
    erdos_renyi_edges(n, p; rng=Random.default_rng())

Generate an Erdős-Rényi random graph G(n, p) and return a vector of
edge tuples (i, j) where 0 ≤ i < j < n (0-indexed vertices).
"""
function erdos_renyi_edges(n::Int, p::Float64; rng=Random.default_rng())
    edges = Tuple{Int,Int}[]
    for i in 0:(n-2)
        for j in (i+1):(n-1)
            if rand(rng) < p
                push!(edges, (i, j))
            end
        end
    end
    return edges
end

"""
    maxcut_optimal(costs)

Find the maximum cut value by brute force over all bitstrings.
"""
function maxcut_optimal(costs::Vector{Float64})
    return maximum(costs)
end

"""
    generate_er_instance(n, p_edge; rng=Random.default_rng())

Generate one Erdős-Rényi graph and return (edges, costs, num_edges).
"""
function generate_er_instance(n::Int, p_edge::Float64; rng=Random.default_rng())
    edges = erdos_renyi_edges(n, p_edge; rng)
    m = length(edges)
    costs = maxcut_costs(n, edges)
    return (edges=edges, costs=costs, num_edges=m, num_vertices=n)
end


# linear_ramp and linear_ramp_matrix are provided by JuliaQAOA module.


#==============================================================================#
#                          Plotting Utilities                                   #
#==============================================================================#

# Standard figure size for paper-style plots
const FIGURE_SIZE = (800, 600)
const HALF_FIGURE_SIZE = (400, 300)

"""Save figure with timestamp in filename."""
function save_figure(fig, name::String; dir=joinpath(@__DIR__, "output"))
    mkpath(dir)
    save(joinpath(dir, name), fig)
    println("Saved: $(joinpath(dir, name))")
end
