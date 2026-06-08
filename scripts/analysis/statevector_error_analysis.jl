# Quantify the numerical agreement between the Julia `qaoa_statevector`
# simulation and QOKit's reference Python FUR backend.
#
# This mirrors the comparison in test/test_qaoa_simulation.jl but, instead of
# a pass/fail tolerance, it reports the distribution of errors (mean, median,
# max) over many random instances. Run from the repo root:
#
#     julia --project scripts/analysis/statevector_error_analysis.jl
#
# Requires the python/qokitvenv virtualenv with QOKit installed (find_python.jl
# auto-detects it).

include(joinpath(@__DIR__, "..", "..", "find_python.jl"))

using JuliaQAOA
using PythonCall
using Printf
using Statistics
import Random: MersenneTwister

const np = pyimport("numpy")
const nx = pyimport("networkx")
const mc = pyimport("qokit.maxcut")
const qaoa_sim = pyimport("grips.QAOA_simulator")

function make_nx_graph(n::Int, edges::Vector{Tuple{Int,Int}})
    G = nx.Graph()
    G.add_nodes_from(pylist(0:(n-1)))
    for (i, j) in edges
        G.add_edge(i, j)
    end
    return G
end

"""
Compare Julia and Python statevectors for one random G(n, 0.5) MaxCut instance
at QAOA depth `p`. Returns a NamedTuple of error summaries for this instance.
"""
function compare_instance(n::Int, p::Int, seed::Int)
    rng = MersenneTwister(seed)
    edges = Tuple{Int,Int}[]
    for i in 0:(n-2), j in (i+1):(n-1)
        if rand(rng) < 0.5
            push!(edges, (i, j))
        end
    end

    param_rng = MersenneTwister(seed + 100)
    γs = rand(param_rng, p) .* 2π
    βs = rand(param_rng, p) .* 2π

    jl_costs = maxcut_costs(n, edges)
    jl_state = qaoa_statevector(jl_costs, n, γs, βs)

    G = make_nx_graph(n, edges)
    terms = mc.get_maxcut_terms(G)
    py_state = pyconvert(
        Vector{ComplexF64},
        np.array(qaoa_sim.get_state(n, terms, np.array(γs), np.array(βs), simulator_name="python")),
    )

    abs_err = abs.(jl_state .- py_state)            # per-amplitude absolute error
    denom = max.(abs.(py_state), eps())              # guard tiny amplitudes
    rel_err = abs_err ./ denom                        # per-amplitude relative error
    return (
        seed = seed,
        n = n,
        p = p,
        nedges = length(edges),
        norm_diff = sqrt(sum(abs2, jl_state .- py_state)),  # ‖Δ‖₂ over the state
        mean_abs = mean(abs_err),
        median_abs = median(abs_err),
        max_abs = maximum(abs_err),
        mean_rel = mean(rel_err),
        median_rel = median(rel_err),
        max_rel = maximum(rel_err),
    )
end

# Match the test (n=8, p=3) but sweep more seeds, plus a couple of larger sizes.
const CONFIGS = [(8, 3, 0:49), (10, 5, 0:19), (12, 8, 0:9)]

function main()
    rows = NamedTuple[]
    for (n, p, seeds) in CONFIGS
        for seed in seeds
            push!(rows, compare_instance(n, p, seed))
        end
    end

    # Pooled per-amplitude statistics weight every instance equally via its own
    # summary; aggregate the per-instance summaries across the whole sweep.
    println("Per-amplitude absolute error |Q_jl - Q_py|:")
    @printf("  mean   = %.3e\n", mean(r.mean_abs for r in rows))
    @printf("  median = %.3e\n", median(r.median_abs for r in rows))
    @printf("  max    = %.3e\n", maximum(r.max_abs for r in rows))
    println()
    println("Per-amplitude relative error |Q_jl - Q_py| / |Q_py|:")
    @printf("  mean   = %.3e\n", mean(r.mean_rel for r in rows))
    @printf("  median = %.3e\n", median(r.median_rel for r in rows))
    @printf("  max    = %.3e\n", maximum(r.max_rel for r in rows))
    println()
    println("Whole-state 2-norm difference ‖Q_jl - Q_py‖₂:")
    @printf("  mean   = %.3e\n", mean(r.norm_diff for r in rows))
    @printf("  median = %.3e\n", median(r.norm_diff for r in rows))
    @printf("  max    = %.3e\n", maximum(r.norm_diff for r in rows))
    println()

    println("Per-configuration max |abs err| and max ‖Δ‖₂:")
    for (n, p, seeds) in CONFIGS
        sub = filter(r -> r.n == n && r.p == p, rows)
        @printf(
            "  n=%2d p=%d (%d instances): max_abs=%.3e  max_norm=%.3e\n",
            n, p, length(sub),
            maximum(r.max_abs for r in sub),
            maximum(r.norm_diff for r in sub),
        )
    end
    return rows
end

main()
