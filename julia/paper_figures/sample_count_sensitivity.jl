#=
sample_count_sensitivity.jl — Test how S (samples per cost class) affects
SampN+EmpP QAOA performance at large n.

Currently S=20 is used everywhere. Is this optimal? Does S=50 or S=100
improve results? This characterizes the accuracy-cost tradeoff for the paper.

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VALUES = [16, 18, 20]
const S_VALUES = [5, 10, 20, 50, 100]
const NUM_HOMODIST_INSTANCES = 10
const NUM_EVAL_INSTANCES = 5

const ER_P_EDGE = 0.5
const P_VALUES = [1, 3]
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const SEED = 42

#==============================================================================#
#                    HELPERS                                                    #
#==============================================================================#

function compute_empirical_P(instances, max_edges)
    P = zeros(max_edges + 1)
    n = instances[1].num_vertices; num_bs = 1 << n
    for inst in instances
        for x in 1:num_bs
            c = Int(inst.costs[x])
            if c <= max_edges; P[c + 1] += 1.0 / (num_bs * length(instances)); end
        end
    end
    return P
end

function optimize_via_proxy(homodist, P_vals, n, p)
    if p == 1
        γ_range = range(0.02, 2.0, length=GRID_SIZE_P1)
        β_range = range(0.01, π/2 - 0.01, length=GRID_SIZE_P1)
        K = GRID_SIZE_P1^2
        γ_matrix = zeros(K, 1); β_matrix = zeros(K, 1)
        idx = 0
        for γ in γ_range, β in β_range
            idx += 1; γ_matrix[idx, 1] = γ / π; β_matrix[idx, 1] = β / π
        end
        Qs = QAOA_proxy_multi(homodist, γ_matrix, β_matrix; pi_units=true)
        exps = vec(expectation(Qs[end], P_vals, n))
        best_idx = argmax(exps)
        return [γ_matrix[best_idx, 1] * π], [β_matrix[best_idx, 1] * π], exps[best_idx]
    else
        gs = GRID_SIZE_RAMP; K = gs^4
        γ_matrix = zeros(K, p); β_matrix = zeros(K, p)
        idx = 0
        for γ₁ in range(0.02, 0.40, length=gs),
            γ_f in range(0.10, 0.70, length=gs),
            β₁ in range(0.05, 0.45, length=gs),
            β_f in range(0.01, 0.25, length=gs)
            idx += 1
            γs, βs = linear_ramp(γ₁, γ_f, β₁, β_f, p)
            γ_matrix[idx, :] .= γs; β_matrix[idx, :] .= βs
        end
        Qs = QAOA_proxy_multi(homodist, γ_matrix, β_matrix; pi_units=true)
        exps = vec(expectation(Qs[end], P_vals, n))
        best_idx = argmax(exps)
        return γ_matrix[best_idx, :] .* π, β_matrix[best_idx, :] .* π, exps[best_idx]
    end
end

#==============================================================================#
#   MAIN EXPERIMENT                                                            #
#==============================================================================#

println("=" ^ 80)
println("Sample Count Sensitivity: S vs QAOA Performance (ER)")
println("S values: $S_VALUES")
println("=" ^ 80)

all_results = Dict{Tuple{Int,Int,Int}, Float64}()  # (n, p, S) -> mean ratio
all_times = Dict{Tuple{Int,Int}, Float64}()          # (n, S) -> homodist time

for n in N_VALUES
    println("\n" * "=" ^ 40)
    @printf("  n = %d\n", n)
    println("=" ^ 40)

    # Generate instances (same for all S values)
    eval_instances = [generate_er_instance(n, ER_P_EDGE;
        rng=MersenneTwister(SEED + n * 100 + i)) for i in 1:NUM_EVAL_INSTANCES]
    eval_optimal = [maxcut_optimal(inst.costs) for inst in eval_instances]

    homodist_pool = [generate_er_instance(n, ER_P_EDGE;
        rng=MersenneTwister(SEED + n * 200 + i)) for i in 1:NUM_HOMODIST_INSTANCES]
    max_edges = maximum(inst.num_edges for inst in homodist_pool)
    P_emp = compute_empirical_P(homodist_pool, max_edges)

    for S in S_VALUES
        println("\n  S=$S:")

        # Compute sampled homodist
        t_sample = @elapsed begin
            sampled_homodists = map(homodist_pool) do inst
                get_homogeneous_distribution_from_costs_sampled(
                    inst.costs, inst.num_edges, inst.num_vertices, S;
                    max_num_edges=max_edges,
                    rng=MersenneTwister(SEED + hash(inst.num_edges) + n + S * 7))
            end
            sampled_avg = average_distributions(sampled_homodists)
        end
        all_times[(n, S)] = t_sample

        m_hd = size(sampled_avg, 1) - 1
        P_vals = P_emp[1:min(m_hd+1, length(P_emp))]
        if length(P_vals) < m_hd + 1
            P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
        end

        @printf("    Homodist computed in %.2fs\n", t_sample)

        for p in P_VALUES
            best_γ, best_β, _ = optimize_via_proxy(sampled_avg, P_vals, n, p)

            ratios = map(enumerate(eval_instances)) do (i, inst)
                qaoa_expectation(inst.costs, n, best_γ, best_β) / eval_optimal[i]
            end
            all_results[(n, p, S)] = mean(ratios)
            @printf("    p=%d: mean=%.4f  std=%.4f\n", p, mean(ratios), std(ratios))
        end
    end
end

#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY — Mean Approx Ratio by Sample Count (ER)")
println("=" ^ 80)

for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-6s", "n")
    for S in S_VALUES; @printf("  S=%-6d", S); end
    println()
    println("  " * "-" ^ (6 + 10 * length(S_VALUES)))
    for n in N_VALUES
        @printf("  %-6d", n)
        for S in S_VALUES
            @printf("  %-8.4f", all_results[(n, p, S)])
        end
        println()
    end
end

println("\n\nTiming — Homodist Computation (seconds, $NUM_HOMODIST_INSTANCES instances):")
@printf("  %-6s", "n")
for S in S_VALUES; @printf("  S=%-6d", S); end
println()
println("  " * "-" ^ (6 + 10 * length(S_VALUES)))
for n in N_VALUES
    @printf("  %-6d", n)
    for S in S_VALUES
        @printf("  %-8.2f", all_times[(n, S)])
    end
    println()
end

println("\nDone!")
