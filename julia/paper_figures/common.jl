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
