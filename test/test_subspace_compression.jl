#=
Tests for the cost-class projection primitives and the exactly-compressed
QAOA trajectory (subspace_compression.jl).

The central identity under test: one step of the homogeneous proxy with the
same-instance empirical N(c'; d, c) equals "apply one true QAOA layer, then
replace each amplitude by its cost-class mean."
=#

using JuliaQAOA, Test
using Random: MersenneTwister
using LinearAlgebra: norm

# Small ER-style random graph helper (deterministic via seed)
function random_edges(rng, n, prob)
    edges = Tuple{Int,Int}[]
    for i in 0:(n-2), j in (i+1):(n-1)
        rand(rng) < prob && push!(edges, (i, j))
    end
    return edges
end


@testset "project_onto_cost_classes basics" begin
    # P₃ path graph: costs [0, 1, 2, 1, 1, 2, 1, 0], classes M = [2, 4, 2]
    costs = maxcut_costs(3, [(0,1), (1,2)])

    # A state constant on classes lies in the subspace: zero residual
    Qin = ComplexF64[0.1 + 0.2im, -0.3im, 0.5]
    state = reconstruct_from_cost_classes(Qin, costs)
    proj = project_onto_cost_classes(state, costs)
    @test proj.counts == [2, 4, 2]
    @test proj.Q ≈ Qin
    @test proj.residual_norm ≈ 0 atol=1e-15

    # Generic state: Pythagoras ‖state‖² = ‖P state‖² + ‖(I−P) state‖²
    rng = MersenneTwister(7)
    state = randn(rng, ComplexF64, 8)
    proj = project_onto_cost_classes(state, costs)
    projected_norm2 = sum(proj.counts .* abs2.(proj.Q))
    @test projected_norm2 + proj.residual_norm^2 ≈ norm(state)^2

    # Idempotence: projecting the reconstruction changes nothing
    reproj = project_onto_cost_classes(reconstruct_from_cost_classes(proj.Q, costs), costs)
    @test reproj.Q ≈ proj.Q
    @test reproj.residual_norm ≈ 0 atol=1e-14

    # num_costs padding for unattained costs (e.g. triangle graph never cuts 1 or 3 edges)
    tri_costs = maxcut_costs(3, [(0,1), (1,2), (0,2)])
    proj = project_onto_cost_classes(fill(0.5 + 0.0im, 8), tri_costs; num_costs=4)
    @test proj.counts == [2, 0, 6, 0]
    @test proj.Q[2] == 0 && proj.Q[4] == 0
end


@testset "compressed trajectory = proxy with empirical N (Theorem 1)" begin
    rng = MersenneTwister(42)
    for n in (6, 8)
        edges = random_edges(rng, n, 0.5)
        m = length(edges)
        costs = maxcut_costs(n, edges)

        γs = [0.3, -0.7, 1.1]
        βs = [0.2, 0.5, -0.4]

        traj = compressed_qaoa_trajectory(costs, n, γs, βs; num_costs=m + 1)
        N = get_homogeneous_distribution_from_costs_direct(costs, m, n)
        Qs_proxy = QAOA_proxy_basic(N, γs, βs)

        # The iterates agree on every attained cost class. Unattained classes
        # are inert: the proxy seeds them with the uniform amplitude at layer 0
        # (where the projection gives 0), but the empirical N has zero
        # rows/columns there, so both evolutions hold them at 0 from layer 1 on.
        attained = traj.counts .> 0
        @test traj.Qs[1][attained] ≈ Qs_proxy[1][attained] rtol=1e-12
        for ℓ in 2:(length(γs) + 1)
            @test traj.Qs[ℓ] ≈ Qs_proxy[ℓ] rtol=1e-12
            @test all(iszero, Qs_proxy[ℓ][.!attained])
        end
    end
end


@testset "leakage diagnostics" begin
    rng = MersenneTwister(3)
    n = 7
    edges = random_edges(rng, n, 0.5)
    costs = maxcut_costs(n, edges)
    p = 5
    γs = 0.8 .* rand(rng, p)
    βs = 0.6 .* rand(rng, p)

    traj = compressed_qaoa_trajectory(costs, n, γs, βs)

    # The compressed norm decays exactly by the per-layer leakage
    for ℓ in 1:p
        @test traj.compressed_norm[ℓ + 1]^2 + traj.leakage[ℓ]^2 ≈ traj.compressed_norm[ℓ]^2
    end

    # Telescoping bound (Theorem 2): ‖ψ_p − φ_p‖ ≤ Σ λ_ℓ
    for ℓ in 1:p
        @test traj.distance[ℓ + 1] <= sum(traj.leakage[1:ℓ]) + 1e-12
    end

    # β = 0 leaves the phase gate alone, which is diagonal in cost classes:
    # the subspace is exactly invariant and nothing leaks
    traj0 = compressed_qaoa_trajectory(costs, n, [1.3, 0.4], [0.0, 0.0])
    @test all(abs.(traj0.leakage) .< 1e-14)
    @test traj0.overlap[end] ≈ 1 atol=1e-12

    # Layer 0 entries describe the uniform state, which is in the subspace
    @test traj.compressed_norm[1] ≈ 1
    @test traj.overlap[1] ≈ 1
    @test traj.distance[1] ≈ 0 atol=1e-15
end
