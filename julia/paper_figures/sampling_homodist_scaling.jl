#=
sampling_homodist_scaling.jl — Test sampling-based homodist at larger n
where exact computation becomes expensive.

Shows the scalability advantage of sampling: O(S × (m+1) × 2^n) vs O(2^(2n)).

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VALUES = [12, 14, 16, 18]
const SAMPLES_PER_COST = 10
const NUM_INSTANCES = 5
const SEED = 42

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

#==============================================================================#
#                    HELPERS                                                    #
#==============================================================================#

function generate_instances(graph_type::String, n::Int, num::Int; rng=Random.default_rng())
    map(1:num) do _
        if graph_type == "ER"
            generate_er_instance(n, ER_P_EDGE; rng)
        elseif graph_type == "BA"
            generate_ba_instance(n, BA_M_ATTACH; rng)
        elseif graph_type == "WS"
            generate_ws_instance(n, WS_K, WS_P_REWIRE; rng)
        end
    end
end

#==============================================================================#
#   SCALING TEST                                                               #
#==============================================================================#

println("=" ^ 80)
println("Sampling-Based Homodist: Scaling with n")
println("=" ^ 80)
println()

graph_types = ["ER", "BA", "WS"]

for gt in graph_types
    println("--- $gt ---")
    @printf("  %-4s  %-8s  %-12s  %-12s  %-8s  %-10s\n",
        "n", "m(avg)", "t_exact(s)", "t_sampled(s)", "speedup", "MSE")

    for n in N_VALUES
        rng = MersenneTwister(SEED + hash(gt) + n)
        instances = generate_instances(gt, n, NUM_INSTANCES; rng)
        max_edges = maximum(inst.num_edges for inst in instances)

        # Time exact homodist
        if n <= 18  # skip exact for very large n
            t_exact = @elapsed begin
                exact_homodists = map(instances) do inst
                    get_homogeneous_distribution_from_costs_direct(
                        inst.costs, inst.num_edges, inst.num_vertices;
                        max_num_edges=max_edges)
                end
                exact_avg = average_distributions(exact_homodists)
            end
        end

        # Time sampled homodist
        t_sampled = @elapsed begin
            sampled_homodists = map(instances) do inst
                get_homogeneous_distribution_from_costs_sampled(
                    inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
                    max_num_edges=max_edges, rng=MersenneTwister(rand(rng, UInt64)))
            end
            sampled_avg = average_distributions(sampled_homodists)
        end

        # Compute MSE if exact is available
        if n <= 18
            sampled_matched, exact_matched = pad_to_match(sampled_avg, exact_avg)
            mse = sum((sampled_matched .- exact_matched).^2) / length(exact_matched)
            speedup = t_exact / max(t_sampled, 1e-6)
            @printf("  %-4d  %-8.0f  %-12.3f  %-12.3f  %-8.1f  %-10.2e\n",
                n, mean(inst.num_edges for inst in instances), t_exact, t_sampled, speedup, mse)
        else
            @printf("  %-4d  %-8.0f  %-12s  %-12.3f  %-8s  %-10s\n",
                n, mean(inst.num_edges for inst in instances), "skipped", t_sampled, "—", "—")
        end
    end
    println()
end

println("Done!")
