#=
dense_nonER_investigation.jl — Test proxy methods on dense non-ER graphs
where PaperProxy is known to fail (p_eff > 0.5).

PaperProxy catastrophically fails at p_eff>0.5 (parameter-sweep entry).
SampN+EmpP doesn't use the analytical formula, so it might work here.

Graph configs:
  - BA(m=4): dense BA, ~2m edges per node → p_eff ~0.6-0.8 for small n
  - WS(k=6,p=0.5): dense WS with high rewiring → p_eff ~0.5-0.6

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VALUES = [12, 14, 16]
const N_SOURCE = 9
const NUM_HOMODIST_INSTANCES = 20
const NUM_EVAL_INSTANCES = 15
const SAMPLES_PER_COST = 20
const SEED = 42

const P_VALUES = [1, 3, 5]
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10

# Dense graph configurations
const CONFIGS = [
    ("BA(m=4)", "BA", 4, 0, 0.0),       # m_attach=4
    ("WS(k=6,p=0.5)", "WS", 0, 6, 0.5), # k=6, p_rewire=0.5
    ("BA(m=2)", "BA", 2, 0, 0.0),        # baseline sparse BA
    ("WS(k=4,p=0.3)", "WS", 0, 4, 0.3),  # baseline sparse WS
]

#==============================================================================#
#                    HELPERS                                                    #
#==============================================================================#

function generate_config_instances(config, n, num; rng=Random.default_rng())
    _, gtype, m_attach, k, p_rewire = config
    map(1:num) do _
        if gtype == "BA"
            generate_ba_instance(n, m_attach; rng)
        else
            generate_ws_instance(n, k, p_rewire; rng)
        end
    end
end

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
println("Dense Non-ER Graph Investigation")
println("=" ^ 80)

methods = ["Transfer", "PaperProxy", "SampN+EmpP"]

for n_target in N_VALUES
    println("\n" * "=" ^ 50)
    println("  n = $n_target")
    println("=" ^ 50)

    for config in CONFIGS
        label, gtype, m_attach, k, p_rewire = config
        println("\n  --- $label (n=$n_target) ---")

        # Generate instances
        eval_instances = generate_config_instances(config, n_target, NUM_EVAL_INSTANCES;
            rng=MersenneTwister(SEED + hash(config) + n_target))
        eval_optimal = [maxcut_optimal(inst.costs) for inst in eval_instances]
        avg_peff = mean(effective_edge_probability(inst.num_edges, n_target) for inst in eval_instances)
        avg_m = mean(inst.num_edges for inst in eval_instances)
        @printf("  avg m=%.0f, avg p_eff=%.3f\n", avg_m, avg_peff)

        # Source instances for transfer
        source_instances = generate_config_instances(config, N_SOURCE, 10;
            rng=MersenneTwister(SEED + hash(config) + 888))

        # Homodist instances
        hd_instances = generate_config_instances(config, n_target, NUM_HOMODIST_INSTANCES;
            rng=MersenneTwister(SEED + hash(config) + n_target * 200))
        max_edges = maximum(inst.num_edges for inst in hd_instances)

        for p in P_VALUES
            # === TRANSFER ===
            if p == 1
                src_params = map(source_instances) do inst
                    best_exp = -Inf; best_γ = 0.0; best_β = 0.0
                    for γ in range(0.02, 2.0, length=80), β in range(0.01, π/2, length=80)
                        e = qaoa_expectation(inst.costs, N_SOURCE, [γ], [β])
                        if e > best_exp; best_exp = e; best_γ = γ; best_β = β; end
                    end
                    (γ=best_γ, β=best_β)
                end
                med_γ = [median([sp.γ for sp in src_params])]
                med_β = [median([sp.β for sp in src_params])]
            else
                src_ramp = map(source_instances) do inst
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
                med_γ = γs_t .* π; med_β = βs_t .* π
            end
            transfer_ratios = map(enumerate(eval_instances)) do (i, inst)
                qaoa_expectation(inst.costs, n_target, med_γ, med_β) / eval_optimal[i]
            end

            # === PAPERPROXY ===
            pp_ratios = map(enumerate(eval_instances)) do (i, inst)
                p_eff = effective_edge_probability(inst.num_edges, n_target)
                proxy = PaperProxy(inst.num_edges, n_target, p_eff)
                homodist = cpu_compute_homodist(proxy)
                P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, n_target, p)
                qaoa_expectation(inst.costs, n_target, best_γ, best_β) / eval_optimal[i]
            end

            # === SAMPLED HOMODIST ===
            sampled_homodists = map(hd_instances) do inst
                get_homogeneous_distribution_from_costs_sampled(
                    inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
                    max_num_edges=max_edges,
                    rng=MersenneTwister(SEED + hash(config) + hash(inst.num_edges)))
            end
            sampled_avg = average_distributions(sampled_homodists)
            P_emp = compute_empirical_P(hd_instances, max_edges)
            m_hd = size(sampled_avg, 1) - 1
            P_vals_s = P_emp[1:min(m_hd+1, length(P_emp))]
            if length(P_vals_s) < m_hd + 1
                P_vals_s = vcat(P_vals_s, zeros(m_hd + 1 - length(P_vals_s)))
            end
            best_γ_s, best_β_s, _ = optimize_via_proxy(sampled_avg, P_vals_s, n_target, p)
            samp_ratios = map(enumerate(eval_instances)) do (i, inst)
                qaoa_expectation(inst.costs, n_target, best_γ_s, best_β_s) / eval_optimal[i]
            end

            @printf("  p=%d:  Transfer=%.4f  PaperProxy=%.4f  SampN+EmpP=%.4f\n",
                p, mean(transfer_ratios), mean(pp_ratios), mean(samp_ratios))
        end
    end
end

println("\nDone!")
