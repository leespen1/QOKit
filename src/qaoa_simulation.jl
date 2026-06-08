#=
qaoa_simulation.jl — Real QAOA statevector simulation (FUR-style).

Provides full 2^n statevector simulation of QAOA circuits using the
Fast Unitary Rotation approach with qubit-by-qubit X-mixer application.
=#


"""
    maxcut_costs(n, edges)

Compute the MaxCut cost for every bitstring x ∈ {0, …, 2^n - 1}.
Returns a Vector{Float64} of length 2^n.

Cost of a cut = number of edges (i,j) where x[i] ≠ x[j].
Vertices are 0-indexed: edges should be tuples (i, j) with 0 ≤ i < j < n.
"""
function maxcut_costs(n::Int, edges::Vector{Tuple{Int,Int}})
    num_bitstrings = 1 << n
    costs = zeros(Float64, num_bitstrings)
    for x in 0:(num_bitstrings - 1)
        c = 0
        for (i, j) in edges
            if ((x >> i) & 1) != ((x >> j) & 1)
                c += 1
            end
        end
        costs[x + 1] = Float64(c)
    end
    return costs
end


"""
    apply_phase_gate!(state, costs, γ)

Apply the phase gate e^{-iγC/2} to the statevector in place.
Each amplitude state[x+1] is multiplied by exp(-iγ·c(x)/2).

Convention: matches QOKit's `exp(-0.5j * gamma * hc_diag)`.
"""
function apply_phase_gate!(state::Vector{Complex{T}}, costs::Vector{T}, γ::Real) where T <: AbstractFloat
    @inbounds for x in eachindex(state)
        state[x] *= cis(T(-γ * costs[x] / 2))
    end
end


"""
    apply_x_mixer!(state, β, n)

Apply the X-mixer e^{-iβΣ_j X_j} using the FUR approach.
For each qubit j, pair up amplitudes that differ only in bit j and apply:
    [cos(β)    -i·sin(β)]
    [-i·sin(β)   cos(β) ]

This is applied qubit-by-qubit, equivalent to ∏_j e^{-iβX_j}.
Convention: matches QOKit's `furx_all(sv, beta, n_qubits)`.
"""
function apply_x_mixer!(state::Vector{Complex{T}}, β::Real, n::Int) where T <: AbstractFloat
    cosβ = T(cos(β))
    sinβ = T(sin(β))
    @inbounds for j in 0:(n-1)
        mask = 1 << j
        for x in 0:(length(state) - 1)
            if (x & mask) == 0  # Only process each pair once
                y = x | mask    # y differs from x at bit j
                ax = state[x + 1]
                ay = state[y + 1]
                state[x + 1] = cosβ * ax - im * sinβ * ay
                state[y + 1] = -im * sinβ * ax + cosβ * ay
            end
        end
    end
end


"""
    qaoa_statevector(costs, n, γs, βs; return_intermediates=false)

Run full QAOA simulation and return the final statevector.
If return_intermediates=true, return a vector of statevectors after each layer
(with the initial uniform state as the first element).

Precision (Float32/Float64) is inferred from the element type of `costs`.
γs and βs are vectors of length p in raw radians (NOT units of π).
"""
function qaoa_statevector(
    costs::Vector{T}, n::Int,
    γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real};
    return_intermediates::Bool=false
) where T <: AbstractFloat
    @assert length(γs) == length(βs) "γs and βs must have the same length"
    p = length(γs)
    num_bitstrings = 1 << n
    CT = Complex{T}

    # Initial state: uniform superposition
    state = fill(CT(1 / sqrt(T(num_bitstrings))), num_bitstrings)

    states = return_intermediates ? [copy(state)] : Vector{CT}[]

    for ℓ in 1:p
        apply_phase_gate!(state, costs, γs[ℓ])
        apply_x_mixer!(state, βs[ℓ], n)
        if return_intermediates
            push!(states, copy(state))
        end
    end

    return return_intermediates ? states : state
end


"""
    qaoa_expectation(costs, n, γs, βs)

Compute ⟨C⟩ = ⟨ψ(γ,β)|C|ψ(γ,β)⟩ from full QAOA simulation.
"""
function qaoa_expectation(
    costs::Vector{T}, n::Int,
    γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real}
) where T <: AbstractFloat
    state = qaoa_statevector(costs, n, γs, βs)
    return real(sum(abs2.(state) .* costs))
end
