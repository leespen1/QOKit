module JuliaQAOAKernelAbstractionsExt

using KernelAbstractions: KernelAbstractions, @kernel, @index, @Const,
    @localmem, @synchronize, get_backend, synchronize
using JuliaQAOA


# ─── MaxCut cost kernel ──────────────────────────────────────────────────────

@kernel function _maxcut_costs_kernel!(costs, @Const(edge_i), @Const(edge_j),
                                       num_edges)
    idx = @index(Global)
    x = idx - 1  # 0-based bitstring
    c = zero(eltype(costs))
    for e in 1:num_edges
        @inbounds if ((x >> edge_i[e]) & 1) != ((x >> edge_j[e]) & 1)
            c += one(eltype(costs))
        end
    end
    @inbounds costs[idx] = c
end


# ─── Phase separator kernel ─────────────────────────────────────────────────

@kernel function _phase_gate_kernel!(state, @Const(costs), γ_half)
    i = @index(Global)
    @inbounds state[i] *= cis(-γ_half * costs[i])
end


"""
    gpu_apply_phase_gate!(state, costs, γ)

Apply the phase gate e^{-iγC/2} to the GPU statevector in place.
Convention: matches QOKit's `exp(-0.5j * gamma * hc_diag)`.
"""
function JuliaQAOA.gpu_apply_phase_gate!(state::AbstractVector{<:Complex},
                               costs::AbstractVector{<:Real}, γ::Real)
    backend = get_backend(state)
    T = real(eltype(state))
    _phase_gate_kernel!(backend)(state, costs, T(γ / 2), ndrange=length(state))
    synchronize(backend)
end


# ─── X-mixer kernel (single qubit) ──────────────────────────────────────────

@kernel function _x_mixer_qubit_kernel!(state, cosβ, sinβ, bit_mask)
    tid = @index(Global)
    t = tid - 1  # 0-based thread index
    x = (t & (bit_mask - 1)) | ((t & ~(bit_mask - 1)) << 1)
    y = x | bit_mask

    @inbounds begin
        ax = state[x + 1]
        ay = state[y + 1]
        state[x + 1] = cosβ * ax - im * sinβ * ay
        state[y + 1] = -im * sinβ * ax + cosβ * ay
    end
end


"""
    gpu_apply_x_mixer!(state, β, n)

Apply the X-mixer e^{-iβΣ_j X_j} using the FUR approach on GPU.
One kernel launch per qubit; within each launch, all 2^(n-1) pairs
are processed in parallel.

Convention: matches QOKit's `furx_all(sv, beta, n_qubits)`.
"""
function JuliaQAOA.gpu_apply_x_mixer!(state::AbstractVector{<:Complex}, β::Real, n::Int)
    backend = get_backend(state)
    T = real(eltype(state))
    cosβ = T(cos(β))
    sinβ = T(sin(β))
    n_pairs = length(state) ÷ 2

    for j in 0:(n-1)
        bit_mask = 1 << j
        _x_mixer_qubit_kernel!(backend)(state, cosβ, sinβ, bit_mask,
                                        ndrange=n_pairs)
        synchronize(backend)
    end
end


# ─── Expectation kernel ─────────────────────────────────────────────────────

@kernel function _weighted_prob_kernel!(out, @Const(state), @Const(costs))
    i = @index(Global)
    @inbounds out[i] = real(conj(state[i]) * state[i]) * costs[i]
end


# ─── Public API ──────────────────────────────────────────────────────────────

"""
    gpu_qaoa_statevector(costs_gpu, n, γs, βs)

Run full QAOA simulation on GPU and return the device statevector.

`costs_gpu` must be a GPU array (e.g. CuArray, oneArray) of length 2^n.
Precision (Float32 or Float64) is inferred from the element type of `costs_gpu`.
γs and βs are vectors of length p in raw radians (NOT units of π).
"""
function JuliaQAOA.gpu_qaoa_statevector(
    costs_gpu::AbstractVector{<:Real}, n::Int,
    γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real}
)
    @assert length(γs) == length(βs) "γs and βs must have the same length"
    p = length(γs)
    num_bitstrings = 1 << n

    backend = get_backend(costs_gpu)
    T = eltype(costs_gpu)
    CT = Complex{T}

    state = KernelAbstractions.allocate(backend, CT, num_bitstrings)
    fill!(state, CT(1 / sqrt(T(num_bitstrings))))

    for ℓ in 1:p
        JuliaQAOA.gpu_apply_phase_gate!(state, costs_gpu, T(γs[ℓ]))
        JuliaQAOA.gpu_apply_x_mixer!(state, T(βs[ℓ]), n)
    end

    return state
end


"""
    gpu_qaoa_expectation(costs_gpu, n, γs, βs)

Compute ⟨C⟩ = ⟨ψ(γ,β)|C|ψ(γ,β)⟩ on GPU.
Returns a scalar (transferred from device).
"""
function JuliaQAOA.gpu_qaoa_expectation(
    costs_gpu::AbstractVector{<:Real}, n::Int,
    γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real}
)
    state = JuliaQAOA.gpu_qaoa_statevector(costs_gpu, n, γs, βs)
    backend = get_backend(state)
    T = real(eltype(state))

    out = KernelAbstractions.allocate(backend, T, length(state))
    _weighted_prob_kernel!(backend)(out, state, costs_gpu, ndrange=length(state))
    synchronize(backend)

    return T(sum(out))
end


"""
    gpu_maxcut_costs(n, edges; backend=nothing, T=Float64)

Compute MaxCut costs for all 2^n bitstrings entirely on GPU.
Use `T=Float32` for backends without Float64 support (Intel, Apple Metal).
If `backend` is not specified, uses KernelAbstractions CPU backend.

Edges are 0-indexed tuples (i, j) with 0 ≤ i < j < n.
"""
function JuliaQAOA.gpu_maxcut_costs(n::Int, edges::Vector{Tuple{Int,Int}};
                          backend=nothing, T::Type{<:AbstractFloat}=Float64)
    if backend === nothing
        backend = KernelAbstractions.CPU()
    end

    num_bitstrings = 1 << n
    num_edges = length(edges)

    ei = Int32[e[1] for e in edges]
    ej = Int32[e[2] for e in edges]

    ei_gpu = KernelAbstractions.allocate(backend, Int32, num_edges)
    ej_gpu = KernelAbstractions.allocate(backend, Int32, num_edges)
    copyto!(ei_gpu, ei)
    copyto!(ej_gpu, ej)

    costs_gpu = KernelAbstractions.allocate(backend, T, num_bitstrings)
    _maxcut_costs_kernel!(backend)(costs_gpu, ei_gpu, ej_gpu, num_edges,
                                   ndrange=num_bitstrings)
    synchronize(backend)

    return costs_gpu
end


# ─── Batched shared-memory X-mixer ─────────────────────────────────────────

include("batched_furx_ka.jl")


# ─── Batched QAOA (uses batched mixer, everything else identical) ──────────

"""
    gpu_qaoa_statevector_batched(costs_gpu, n, γs, βs; group_size=10)

Like `gpu_qaoa_statevector` but uses the batched shared-memory X-mixer.
"""
function JuliaQAOA.gpu_qaoa_statevector_batched(
    costs_gpu::AbstractVector{<:Real}, n::Int,
    γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real};
    group_size::Int=10
)
    @assert length(γs) == length(βs) "γs and βs must have the same length"
    p = length(γs)
    num_bitstrings = 1 << n

    backend = get_backend(costs_gpu)
    T = eltype(costs_gpu)
    CT = Complex{T}

    state = KernelAbstractions.allocate(backend, CT, num_bitstrings)
    fill!(state, CT(1 / sqrt(T(num_bitstrings))))

    for ℓ in 1:p
        JuliaQAOA.gpu_apply_phase_gate!(state, costs_gpu, T(γs[ℓ]))
        JuliaQAOA.gpu_apply_x_mixer_batched!(state, T(βs[ℓ]), n; group_size)
    end

    return state
end


"""
    gpu_qaoa_expectation_batched(costs_gpu, n, γs, βs; group_size=10)

Like `gpu_qaoa_expectation` but uses the batched shared-memory X-mixer.
"""
function JuliaQAOA.gpu_qaoa_expectation_batched(
    costs_gpu::AbstractVector{<:Real}, n::Int,
    γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real};
    group_size::Int=10
)
    state = JuliaQAOA.gpu_qaoa_statevector_batched(costs_gpu, n, γs, βs; group_size)
    backend = get_backend(state)
    T = real(eltype(state))

    out = KernelAbstractions.allocate(backend, T, length(state))
    _weighted_prob_kernel!(backend)(out, state, costs_gpu, ndrange=length(state))
    synchronize(backend)

    return T(sum(out))
end


end # module
