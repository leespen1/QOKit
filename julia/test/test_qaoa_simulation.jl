using JuliaQAOA, Test
using PythonCall
import Random: MersenneTwister

# Set up Python environment: import QOKit modules
const np = pyimport("numpy")
const nx = pyimport("networkx")
const mc = pyimport("qokit.maxcut")
const qaoa_sim = pyimport("grips.QAOA_simulator")

"""
    make_nx_graph(n, edges)

Create a NetworkX graph from 0-indexed edge tuples, matching the Julia
graph representation used throughout this codebase.
"""
function make_nx_graph(n::Int, edges::Vector{Tuple{Int,Int}})
    G = nx.Graph()
    G.add_nodes_from(pylist(0:(n-1)))
    for (i, j) in edges
        G.add_edge(i, j)
    end
    return G
end

@testset "maxcut_costs matches QOKit" begin
    for seed in 0:4
        rng = MersenneTwister(seed)
        n = 8
        edges = Tuple{Int,Int}[]
        for i in 0:(n-2), j in (i+1):(n-1)
            if rand(rng) < 0.5
                push!(edges, (i, j))
            end
        end

        # Julia costs
        jl_costs = maxcut_costs(n, edges)

        # Python costs via QOKit
        G = make_nx_graph(n, edges)
        terms = mc.get_maxcut_terms(G)
        from_qokit = qaoa_sim.get_simulator(n, terms, simulator_name="python")
        py_costs = pyconvert(Vector{Float64}, np.array(from_qokit.get_cost_diagonal()))

        @test jl_costs == py_costs
    end
end

@testset "qaoa_statevector matches QOKit" begin
    for seed in 0:4
        rng = MersenneTwister(seed)
        n = 8
        edges = Tuple{Int,Int}[]
        for i in 0:(n-2), j in (i+1):(n-1)
            if rand(rng) < 0.5
                push!(edges, (i, j))
            end
        end

        # Random QAOA parameters (in radians)
        p = 3
        param_rng = MersenneTwister(seed + 100)
        γs = rand(param_rng, p) .* 2π
        βs = rand(param_rng, p) .* 2π

        # Julia statevector
        jl_costs = maxcut_costs(n, edges)
        jl_state = qaoa_statevector(jl_costs, n, γs, βs)

        # Python statevector via QOKit
        G = make_nx_graph(n, edges)
        terms = mc.get_maxcut_terms(G)
        py_gamma = np.array(γs)
        py_beta = np.array(βs)
        py_state = pyconvert(
            Vector{ComplexF64},
            np.array(qaoa_sim.get_state(n, terms, py_gamma, py_beta, simulator_name="python"))
        )

        @test jl_state ≈ py_state atol=1e-18
    end
end

@testset "qaoa_expectation matches QOKit" begin
    for seed in 0:4
        rng = MersenneTwister(seed)
        n = 8
        edges = Tuple{Int,Int}[]
        for i in 0:(n-2), j in (i+1):(n-1)
            if rand(rng) < 0.5
                push!(edges, (i, j))
            end
        end

        p = 3
        param_rng = MersenneTwister(seed + 100)
        γs = rand(param_rng, p) .* 2π
        βs = rand(param_rng, p) .* 2π

        # Julia expectation
        jl_costs = maxcut_costs(n, edges)
        jl_exp = qaoa_expectation(jl_costs, n, γs, βs)

        # Python expectation via QOKit
        G = make_nx_graph(n, edges)
        terms = mc.get_maxcut_terms(G)
        py_gamma = np.array(γs)
        py_beta = np.array(βs)
        py_exp = pyconvert(Float64, qaoa_sim.get_expectation(
            n, terms, py_gamma, py_beta, simulator_name="python"
        ))

        @test jl_exp ≈ py_exp atol=1e-10
    end
end

@testset "qaoa_statevector intermediate states match QOKit layer-by-layer" begin
    rng = MersenneTwister(42)
    n = 6
    edges = Tuple{Int,Int}[]
    for i in 0:(n-2), j in (i+1):(n-1)
        if rand(rng) < 0.5
            push!(edges, (i, j))
        end
    end

    p = 5
    param_rng = MersenneTwister(142)
    γs = rand(param_rng, p) .* 2π
    βs = rand(param_rng, p) .* 2π

    # Julia: get all intermediate states
    jl_costs = maxcut_costs(n, edges)
    jl_states = qaoa_statevector(jl_costs, n, γs, βs; return_intermediates=true)

    # Python: run layer by layer to get intermediates
    G = make_nx_graph(n, edges)
    terms = mc.get_maxcut_terms(G)

    for layer in 1:p
        py_state = pyconvert(
            Vector{ComplexF64},
            np.array(qaoa_sim.get_state(
                n, terms,
                np.array(γs[1:layer]),
                np.array(βs[1:layer]),
                simulator_name="python"
            ))
        )
        # jl_states[1] is initial state, jl_states[layer+1] is state after layer layers
        @test jl_states[layer + 1] ≈ py_state atol=1e-12
    end
end
