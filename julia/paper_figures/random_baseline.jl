#=
random_baseline.jl — Compute random partition baseline approximation ratios
for all experimental configurations.

Random baseline = mean(costs) / max(costs) = E[c(random partition)] / c_opt.
For MaxCut, E[c(random)] = m/2 for unbiased partitions. For biased partitions
with vertex probability q, E[c] = 2q(1-q)m for ER graphs.

This contextualizes all prior QAOA improvements: if the random baseline is
already 0.9, a proxy achieving 0.92 is less impressive than it looks.

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VALUES = [12, 14, 16, 18, 20, 22]
const NUM_INSTANCES = 20
const SEED = 42

# Graph configurations
const ER_P_EDGE = 0.5
const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3

# ER with p != 0.5
const ER_P_EDGES = [0.1, 0.2, 0.3, 0.5, 0.7]

#==============================================================================#
#                          HELPERS                                              #
#==============================================================================#

"""Compute random baseline: mean(costs)/max(costs) for unbiased partition."""
function random_baseline(costs::Vector{Float64})
    return mean(costs) / maximum(costs)
end

"""Compute biased baseline: for vertex probability q, E[c] = Σ_edges 2q(1-q).
For general graphs, compute exactly as average over biased random partitions."""
function biased_baseline_exact(costs::Vector{Float64}, n::Int; num_samples::Int=10000, rng=Random.default_rng())
    c_opt = maximum(costs)
    best_ratio = 0.0
    # Try different bias levels
    for q in range(0.1, 0.9, length=17)
        total = 0.0
        for _ in 1:num_samples
            # Generate biased random bitstring
            x = 0
            for bit in 0:(n-1)
                if rand(rng) < q
                    x |= (1 << bit)
                end
            end
            total += costs[x + 1]
        end
        ratio = (total / num_samples) / c_opt
        best_ratio = max(best_ratio, ratio)
    end
    return best_ratio
end

"""For ER(n, p_edge), the unbiased random baseline is theoretically m/2 / c_opt."""
function theoretical_random_baseline_er(num_edges::Int, c_opt::Float64)
    return (num_edges / 2.0) / c_opt
end

#==============================================================================#
#                          MAIN EXPERIMENT                                      #
#==============================================================================#

println("=" ^ 80)
println("Random Partition Baseline Approximation Ratios")
println("=" ^ 80)

# --- Part 1: Standard graph types (ER(0.5), BA, WS) ---

println("\n--- Part 1: ER(0.5), BA(m=2), WS(k=4, p=0.3) ---\n")
@printf("%-6s %-4s  %-8s %-14s %-14s  %-8s\n",
    "Graph", "n", "m(avg)", "Random(unb)", "Random(biased)", "c_opt(avg)")
println("-" ^ 70)

for graph_type in ["ER", "BA", "WS"]
    for n in N_VALUES
        # Skip large n for non-ER (costs array is O(2^n))
        if n > 20 && graph_type != "ER"
            continue
        end

        rng = MersenneTwister(SEED + hash(graph_type) + n * 100)
        instances = map(1:NUM_INSTANCES) do _
            if graph_type == "ER"
                generate_er_instance(n, ER_P_EDGE; rng)
            elseif graph_type == "BA"
                generate_ba_instance(n, BA_M_ATTACH; rng)
            else
                generate_ws_instance(n, WS_K, WS_P_REWIRE; rng)
            end
        end

        m_avg = mean(inst.num_edges for inst in instances)
        c_opts = [maxcut_optimal(inst.costs) for inst in instances]
        unbiased = [random_baseline(inst.costs) for inst in instances]
        biased = [biased_baseline_exact(inst.costs, n; rng=MersenneTwister(SEED + hash(graph_type) + n * 300 + i))
                  for (i, inst) in enumerate(instances)]

        @printf("%-6s %-4d  %-8.1f %.4f±%.4f  %.4f±%.4f  %-8.1f\n",
            graph_type, n, m_avg,
            mean(unbiased), std(unbiased),
            mean(biased), std(biased),
            mean(c_opts))
    end
end

# --- Part 2: ER with different p_edge values ---

println("\n\n--- Part 2: ER(n, p_edge) for varying p_edge ---\n")
@printf("%-6s %-4s  %-8s %-14s %-14s  %-8s  %-8s\n",
    "p_edge", "n", "m(avg)", "Random(unb)", "Random(biased)", "c_opt", "m/2/copt")
println("-" ^ 80)

for p_edge in ER_P_EDGES
    for n in N_VALUES
        # Skip large n for very sparse or dense graphs
        if n > 20 && p_edge != 0.5
            continue
        end
        if n > 18 && p_edge <= 0.1
            continue
        end

        rng = MersenneTwister(SEED + hash(p_edge) + n * 100)
        instances = map(1:NUM_INSTANCES) do _
            generate_er_instance(n, p_edge; rng)
        end

        m_avg = mean(inst.num_edges for inst in instances)
        c_opts = [maxcut_optimal(inst.costs) for inst in instances]
        unbiased = [random_baseline(inst.costs) for inst in instances]
        biased = [biased_baseline_exact(inst.costs, n; rng=MersenneTwister(SEED + hash(p_edge) + n * 300 + i))
                  for (i, inst) in enumerate(instances)]
        theoretical = [theoretical_random_baseline_er(inst.num_edges, maxcut_optimal(inst.costs))
                       for inst in instances]

        @printf("%-6.1f %-4d  %-8.1f %.4f±%.4f  %.4f±%.4f  %-8.1f  %.4f\n",
            p_edge, n, m_avg,
            mean(unbiased), std(unbiased),
            mean(biased), std(biased),
            mean(c_opts),
            mean(theoretical))
    end
end

# --- Part 3: Compare random baseline to prior QAOA results ---

println("\n\n--- Part 3: Context for Prior Results ---\n")
println("Key comparisons (from prior experiments):")
println()
println("At n=22, p=3 on ER(0.5):")
println("  SampN+EmpP:  0.881 (headline result)")
println("  Transfer:    0.844")
println("  Random baseline should be computed above for comparison.")
println()
println("At n=18, p=3 on ER(0.5):")
println("  SampN+EmpP:  0.865")
println("  Transfer:    0.846")
println()
println("The gap between QAOA methods and random baseline is the 'real' improvement.")

println("\nDone!")
