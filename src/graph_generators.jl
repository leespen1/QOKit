#=
graph_generators.jl — Random graph families as MaxCut edge lists.

All generators return a Vector{Tuple{Int,Int}} of edges (i, j) with
0 ≤ i < j < n (0-indexed vertices), the format `maxcut_costs` expects.
Pass an explicit rng for reproducibility; every experiment script should.

Erdős–Rényi keeps the original pair-loop implementation (moved here from
scripts/paper_figures/common.jl) so that existing seeds reproduce the same
instances. The other families wrap Graphs.jl generators.
=#

using Graphs: barabasi_albert, watts_strogatz, random_regular_graph, edges, src, dst
using Random: AbstractRNG, default_rng

export erdos_renyi_edges, barabasi_albert_edges, watts_strogatz_edges, random_regular_edges


# Convert a Graphs.jl graph to sorted 0-indexed edge tuples
function _edge_tuples(g)
    out = [(min(src(e), dst(e)) - 1, max(src(e), dst(e)) - 1) for e in edges(g)]
    return sort!(out)
end


"""
    erdos_renyi_edges(n, p; rng=default_rng())

Generate an Erdős–Rényi G(n, p) graph: each of the n(n-1)/2 possible edges is
included independently with probability `p`.
"""
function erdos_renyi_edges(n::Integer, p::Real; rng::AbstractRNG=default_rng())
    edge_list = Tuple{Int,Int}[]
    for i in 0:(n-2)
        for j in (i+1):(n-1)
            if rand(rng) < p
                push!(edge_list, (i, j))
            end
        end
    end
    return edge_list
end


"""
    barabasi_albert_edges(n, k; rng=default_rng())

Generate a Barabási–Albert preferential-attachment graph: starting from `k`
vertices, each new vertex attaches to `k` existing vertices with probability
proportional to their degree. Produces m = k(n-k) edges.
"""
function barabasi_albert_edges(n::Integer, k::Integer; rng::AbstractRNG=default_rng())
    return _edge_tuples(barabasi_albert(n, k; rng))
end


"""
    watts_strogatz_edges(n, k, β; rng=default_rng())

Generate a Watts–Strogatz small-world graph: a ring lattice where each vertex
connects to its `k` nearest neighbors (`k` even), with each edge rewired to a
random vertex with probability `β`. Produces m = nk/2 edges.
"""
function watts_strogatz_edges(n::Integer, k::Integer, β::Real; rng::AbstractRNG=default_rng())
    return _edge_tuples(watts_strogatz(n, k, β; rng))
end


"""
    random_regular_edges(n, k; rng=default_rng())

Generate a uniformly random `k`-regular graph on `n` vertices (`n*k` even).
Produces m = nk/2 edges.
"""
function random_regular_edges(n::Integer, k::Integer; rng::AbstractRNG=default_rng())
    return _edge_tuples(random_regular_graph(n, k; rng))
end
