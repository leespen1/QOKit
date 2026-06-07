#=
low_peff_large_n.jl — Test PaperProxy with low p_eff (0.2) on non-ER graphs
at n=14-18 and p=1-8.

The paperproxy-advantage investigation found p_eff=0.1-0.2 works best for non-ER
at n=12. Does this generalize to larger n and higher depths?

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VALUES = [12, 14, 16, 18]
const N_SOURCE = 9
const NUM_EVAL_INSTANCES = 10
const SEED = 42

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

const P_VALUES = [1, 3, 5, 8]
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10

const LOW_PEFF = 0.2  # the "regularizer" setting

#==============================================================================#
#                    HELPERS                                                    #
#==============================================================================#

function generate_instances(graph_type::String, n::Int, num::Int; rng=Random.default_rng())
    map(1:num) do _
        if graph_type == "ER"; generate_er_instance(n, ER_P_EDGE; rng)
        elseif graph_type == "BA"; generate_ba_instance(n, BA_M_ATTACH; rng)
        elseif graph_type == "WS"; generate_ws_instance(n, WS_K, WS_P_REWIRE; rng)
        else; error("Unknown: $graph_type"); end
    end
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
println("PaperProxy with Low p_eff on Non-ER Graphs")
println("=" ^ 80)

graph_types = ["BA", "WS"]
methods = ["Transfer", "PP(natural)", "PP(0.2)"]

# Source instances for transfer
source_instances = Dict(gt => generate_instances(gt, N_SOURCE, 10;
    rng=MersenneTwister(SEED + hash(gt) + 888)) for gt in graph_types)

for n_target in N_VALUES
    println("\n" * "=" ^ 40)
    println("  n = $n_target")
    println("=" ^ 40)

    eval_instances = Dict(gt => generate_instances(gt, n_target, NUM_EVAL_INSTANCES;
        rng=MersenneTwister(SEED + hash(gt) + n_target * 100)) for gt in graph_types)
    eval_optimal = Dict(gt => [maxcut_optimal(inst.costs) for inst in eval_instances[gt]]
        for gt in graph_types)

    for gt in graph_types
        println("\n--- $gt (n=$n_target) ---")

        for p in P_VALUES
            # Skip p=5,8 with full ramp grid at large n — too slow
            if p > 3 && n_target > 16
                println("  p=$p: (skipping, QAOA eval too slow)")
                continue
            end

            # === TRANSFER ===
            if p == 1
                src_params = map(source_instances[gt]) do inst
                    best_exp = -Inf; best_γ = 0.0; best_β = 0.0
                    for γ in range(0.02, 2.0, length=80), β in range(0.01, π/2, length=80)
                        e = qaoa_expectation(inst.costs, N_SOURCE, [γ], [β])
                        if e > best_exp; best_exp = e; best_γ = γ; best_β = β; end
                    end
                    (γ=best_γ, β=best_β)
                end
                med_γ = [median([sp.γ for sp in src_params])]
                med_β = [median([sp.β for sp in src_params])]
                transfer_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                    qaoa_expectation(inst.costs, n_target, med_γ, med_β) / eval_optimal[gt][i]
                end
            else
                src_ramp = map(source_instances[gt]) do inst
                    best_params = (0.0, 0.0, 0.0, 0.0); best_exp = -Inf; gs = GRID_SIZE_RAMP
                    for γ₁ in range(0.02, 0.40, length=gs),
                        γ_f in range(0.10, 0.70, length=gs),
                        β₁ in range(0.05, 0.45, length=gs),
                        β_f in range(0.01, 0.25, length=gs)
                        γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, p)
                        exp_val = qaoa_expectation(inst.costs, N_SOURCE, γs_pi .* π, βs_pi .* π)
                        if exp_val > best_exp; best_exp = exp_val; best_params = (γ₁, γ_f, β₁, β_f); end
                    end
                    best_params
                end
                med_ramp = Tuple(median([sp[i] for sp in src_ramp]) for i in 1:4)
                γs_t, βs_t = linear_ramp(med_ramp..., p)
                transfer_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                    qaoa_expectation(inst.costs, n_target, γs_t .* π, βs_t .* π) / eval_optimal[gt][i]
                end
            end

            # === PAPERPROXY (natural p_eff) ===
            pp_nat_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                p_eff = effective_edge_probability(inst.num_edges, n_target)
                proxy = PaperProxy(inst.num_edges, n_target, p_eff)
                homodist = cpu_compute_homodist(proxy)
                P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, n_target, p)
                qaoa_expectation(inst.costs, n_target, best_γ, best_β) / eval_optimal[gt][i]
            end

            # === PAPERPROXY (low p_eff = 0.2) ===
            pp_low_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                proxy = PaperProxy(inst.num_edges, n_target, LOW_PEFF)
                homodist = cpu_compute_homodist(proxy)
                P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, n_target, p)
                qaoa_expectation(inst.costs, n_target, best_γ, best_β) / eval_optimal[gt][i]
            end

            @printf("  p=%d:  Transfer=%.4f  PP(nat)=%.4f  PP(0.2)=%.4f\n",
                p, mean(transfer_ratios), mean(pp_nat_ratios), mean(pp_low_ratios))
        end
    end
end

println("\nDone!")
