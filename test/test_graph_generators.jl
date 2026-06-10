#=
Tests for the random-graph edge-list generators (graph_generators.jl).
All generators must return 0-indexed (i, j) tuples with i < j, the format
maxcut_costs expects, and be deterministic given an rng.
=#

using JuliaQAOA, Test
using Random: MersenneTwister

valid_edges(edge_list, n) = all(e -> 0 <= e[1] < e[2] < n, edge_list) &&
                            allunique(edge_list)

degrees(edge_list, n) = begin
    deg = zeros(Int, n)
    for (i, j) in edge_list
        deg[i + 1] += 1
        deg[j + 1] += 1
    end
    deg
end


@testset "edge-list invariants and counts" begin
    n = 14
    rng = MersenneTwister(2)

    er = erdos_renyi_edges(n, 0.5; rng)
    @test valid_edges(er, n)
    @test length(er) <= n * (n - 1) ÷ 2
    @test erdos_renyi_edges(n, 1.0; rng) |> length == n * (n - 1) ÷ 2
    @test erdos_renyi_edges(n, 0.0; rng) |> isempty

    k = 3
    ba = barabasi_albert_edges(n, k; rng)
    @test valid_edges(ba, n)
    @test length(ba) == k * (n - k)

    ws = watts_strogatz_edges(n, 4, 0.2; rng)
    @test valid_edges(ws, n)
    @test length(ws) == n * 4 ÷ 2
    # β = 0 is the unrewired ring lattice: every vertex has degree k
    ws0 = watts_strogatz_edges(n, 4, 0.0; rng)
    @test all(==(4), degrees(ws0, n))

    rr = random_regular_edges(n, k; rng)
    @test valid_edges(rr, n)
    @test all(==(k), degrees(rr, n))
end


@testset "determinism given a seed" begin
    for gen in (
        rng -> erdos_renyi_edges(10, 0.4; rng),
        rng -> barabasi_albert_edges(10, 2; rng),
        rng -> watts_strogatz_edges(10, 4, 0.3; rng),
        rng -> random_regular_edges(10, 3; rng),
    )
        @test gen(MersenneTwister(11)) == gen(MersenneTwister(11))
    end
end


@testset "compatibility with maxcut_costs" begin
    rng = MersenneTwister(5)
    n = 8
    for edge_list in (
        erdos_renyi_edges(n, 0.5; rng),
        barabasi_albert_edges(n, 2; rng),
        watts_strogatz_edges(n, 4, 0.2; rng),
        random_regular_edges(n, 3; rng),
    )
        costs = maxcut_costs(n, edge_list)
        @test length(costs) == 2^n
        @test all(0 .<= costs .<= length(edge_list))
        # Complement symmetry of MaxCut: c(x) = c(~x)
        @test costs == reverse(costs)
    end
end
