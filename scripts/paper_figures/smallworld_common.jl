#=
smallworld_common.jl — Graph generators for Barabási-Albert and Watts-Strogatz
graphs, plus shared utilities for the small-world investigation.

These are pure-Julia implementations (no NetworkX dependency) that produce
edge lists in the same (i, j) 0-indexed format as erdos_renyi_edges() in common.jl.
=#

include("common.jl")
using Statistics: mean, std, median

#==============================================================================#
#                    Barabási-Albert Preferential Attachment                    #
#==============================================================================#

"""
    barabasi_albert_edges(n, m_attach; rng=Random.default_rng())

Generate a Barabási-Albert preferential attachment graph.

- `n`: number of nodes
- `m_attach`: number of edges each new node attaches (must satisfy 1 ≤ m_attach < n)

Returns a vector of edge tuples (i, j) with 0 ≤ i < j < n.
"""
function barabasi_albert_edges(n::Int, m_attach::Int; rng=Random.default_rng())
    @assert 1 <= m_attach < n "m_attach must satisfy 1 ≤ m_attach < n"

    # Start with a complete graph on (m_attach + 1) nodes
    edges = Set{Tuple{Int,Int}}()
    for i in 0:m_attach
        for j in (i+1):m_attach
            push!(edges, (i, j))
        end
    end

    # Degree array for preferential attachment (0-indexed nodes)
    degree = zeros(Int, n)
    for (i, j) in edges
        degree[i+1] += 1
        degree[j+1] += 1
    end

    # Add nodes one by one
    for new_node in (m_attach + 1):(n - 1)
        # Select m_attach distinct targets with probability proportional to degree
        targets = Set{Int}()
        total_degree = sum(degree[1:new_node])  # only existing nodes
        while length(targets) < m_attach
            r = rand(rng) * total_degree
            cumulative = 0.0
            for node in 0:(new_node - 1)
                cumulative += degree[node + 1]
                if cumulative >= r
                    push!(targets, node)
                    break
                end
            end
        end

        for target in targets
            push!(edges, (min(target, new_node), max(target, new_node)))
            degree[new_node + 1] += 1
            degree[target + 1] += 1
        end
    end

    return sort(collect(edges))
end


#==============================================================================#
#                    Watts-Strogatz Small-World Model                          #
#==============================================================================#

"""
    watts_strogatz_edges(n, k, p_rewire; rng=Random.default_rng())

Generate a Watts-Strogatz small-world graph.

- `n`: number of nodes (arranged in a ring)
- `k`: each node connected to k nearest neighbors (must be even, ≥ 2)
- `p_rewire`: rewiring probability (0 = regular lattice, 1 = random)

Returns a vector of edge tuples (i, j) with 0 ≤ i < j < n.
"""
function watts_strogatz_edges(n::Int, k::Int, p_rewire::Float64; rng=Random.default_rng())
    @assert k >= 2 && iseven(k) "k must be even and ≥ 2"
    @assert k < n "k must be < n"

    # Start with ring lattice: each node i connected to i+1, i+2, ..., i+k/2 (mod n)
    edges = Set{Tuple{Int,Int}}()
    for i in 0:(n-1)
        for offset in 1:(k ÷ 2)
            j = mod(i + offset, n)
            push!(edges, (min(i, j), max(i, j)))
        end
    end

    # Rewire each edge with probability p_rewire
    edges_to_add = Tuple{Int,Int}[]
    edges_to_remove = Tuple{Int,Int}[]

    for i in 0:(n-1)
        for offset in 1:(k ÷ 2)
            j = mod(i + offset, n)
            edge = (min(i, j), max(i, j))
            if rand(rng) < p_rewire && edge in edges
                # Rewire: keep i, choose new target
                for _ in 1:100  # max attempts
                    new_j = rand(rng, 0:(n-1))
                    if new_j != i
                        new_edge = (min(i, new_j), max(i, new_j))
                        if new_edge ∉ edges && new_edge ∉ edges_to_add
                            push!(edges_to_remove, edge)
                            push!(edges_to_add, new_edge)
                            break
                        end
                    end
                end
            end
        end
    end

    for e in edges_to_remove
        delete!(edges, e)
    end
    for e in edges_to_add
        push!(edges, e)
    end

    return sort(collect(edges))
end


#==============================================================================#
#                    Unified Graph Instance Generation                         #
#==============================================================================#

"""
    generate_ba_instance(n, m_attach; rng)

Generate a Barabási-Albert graph and return (edges, costs, num_edges, num_vertices).
"""
function generate_ba_instance(n::Int, m_attach::Int; rng=Random.default_rng())
    edges = barabasi_albert_edges(n, m_attach; rng)
    m = length(edges)
    costs = maxcut_costs(n, edges)
    return (edges=edges, costs=costs, num_edges=m, num_vertices=n)
end

"""
    generate_ws_instance(n, k, p_rewire; rng)

Generate a Watts-Strogatz graph and return (edges, costs, num_edges, num_vertices).
"""
function generate_ws_instance(n::Int, k::Int, p_rewire::Float64; rng=Random.default_rng())
    edges = watts_strogatz_edges(n, k, p_rewire; rng)
    m = length(edges)
    costs = maxcut_costs(n, edges)
    return (edges=edges, costs=costs, num_edges=m, num_vertices=n)
end

"""
    effective_edge_probability(num_edges, num_vertices)

Compute p_eff = m / (n choose 2), the density of the graph.
This is the natural analog of the ER edge probability for non-ER graphs.
"""
function effective_edge_probability(num_edges::Int, num_vertices::Int)
    max_possible = num_vertices * (num_vertices - 1) ÷ 2
    return num_edges / max_possible
end
