#=
common.jl — Shared utilities for paper figure reproduction scripts.

Provides:
  - Erdős-Rényi graph generation (as edge lists)
  - MaxCut cost computation for all 2^n bitstrings
  - Real QAOA statevector simulation (FUR-style qubit-by-qubit X-mixer)
  - Linear ramp schedule generation
  - JuliaQAOA module loading
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JuliaQAOA
using CairoMakie
using Random
using LinearAlgebra

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
    maxcut_costs(n, edges)

Compute the MaxCut cost for every bitstring x ∈ {0, …, 2^n - 1}.
Returns a Vector{Float64} of length 2^n.

Cost of a cut = number of edges (i,j) where x[i] ≠ x[j].
"""
function maxcut_costs(n::Int, edges::Vector{Tuple{Int,Int}})
    num_bitstrings = 1 << n
    costs = zeros(Float64, num_bitstrings)
    for x in 0:(num_bitstrings - 1)
        c = 0
        for (i, j) in edges
            # Check if bits i and j differ
            if ((x >> i) & 1) != ((x >> j) & 1)
                c += 1
            end
        end
        costs[x + 1] = Float64(c)
    end
    return costs
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


#==============================================================================#
#                    Real QAOA Statevector Simulation                           #
#==============================================================================#

"""
    apply_phase_gate!(state, costs, γ)

Apply the phase gate e^{-iγC} to the statevector in place.
Each amplitude state[x+1] is multiplied by exp(-iγ·c(x)).
"""
function apply_phase_gate!(state::Vector{ComplexF64}, costs::Vector{Float64}, γ::Real)
    @inbounds for x in eachindex(state)
        state[x] *= cis(-γ * costs[x])
    end
end

"""
    apply_x_mixer!(state, β, n)

Apply the X-mixer e^{-iβB} where B = Σ_j X_j, using the FUR approach.
For each qubit j, pair up amplitudes that differ only in bit j and apply:
    [cos(β)    -i·sin(β)]
    [-i·sin(β)   cos(β) ]

This is applied qubit-by-qubit, which is equivalent to the tensor product.
"""
function apply_x_mixer!(state::Vector{ComplexF64}, β::Real, n::Int)
    cosβ = cos(β)
    sinβ = sin(β)
    @inbounds for j in 0:(n-1)
        mask = 1 << j
        for x in 0:(length(state) - 1)
            if (x & mask) == 0  # Only process each pair once
                y = x | mask    # y differs from x at bit j
                # Apply 2x2 rotation
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
If return_intermediates=true, return a vector of statevectors after each layer.

γs and βs are vectors of length p (NOT in units of π — raw radians).
"""
function qaoa_statevector(
    costs::Vector{Float64}, n::Int,
    γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real};
    return_intermediates::Bool=false
)
    @assert length(γs) == length(βs) "γs and βs must have the same length"
    p = length(γs)
    num_bitstrings = 1 << n

    # Initial state: uniform superposition
    state = fill(ComplexF64(1.0 / sqrt(num_bitstrings)), num_bitstrings)

    states = return_intermediates ? [copy(state)] : Vector{ComplexF64}[]

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
    costs::Vector{Float64}, n::Int,
    γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real}
)
    state = qaoa_statevector(costs, n, γs, βs)
    return real(sum(abs2.(state) .* costs))
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
