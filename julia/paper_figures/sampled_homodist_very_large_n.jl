#=
sampled_homodist_very_large_n.jl — Test SampN+EmpP at n=20-24 on ER graphs,
where it showed an advantage over Transfer at n>=16.

At these sizes:
  - Exact homodist is INFEASIBLE (O(2^(2n)) = 10^12+ ops)
  - Costs array fits in memory up to n≈25 (2^25 = 33M entries)
  - Only Transfer, PaperProxy, and SampN+EmpP are available

Key question: Does the SampN+EmpP advantage over Transfer continue to grow?

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VALUES = [18, 20, 22]       # target graph sizes (n=24 may be too slow for QAOA eval)
const N_SOURCE = 9                  # source for transfer
const NUM_HOMODIST_INSTANCES = 10   # fewer instances at large n (costs computation is expensive)
const NUM_EVAL_INSTANCES = 5        # fewer eval instances (QAOA simulation is O(2^n))
const SAMPLES_PER_COST = 20

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
println("SampN+EmpP at Very Large n (ER graphs only)")
println("=" ^ 80)

# Source instances for transfer
source_instances = [generate_er_instance(N_SOURCE, ER_P_EDGE;
    rng=MersenneTwister(SEED + 888 + i)) for i in 1:10]

for n_target in N_VALUES
    println("\n" * "=" ^ 40)
    @printf("  n = %d  (2^n = %d, memory ~%.1f MB for costs)\n",
        n_target, 1 << n_target, (1 << n_target) * 8 / 1e6)
    println("=" ^ 40)

    # Generate instances
    println("  Generating $NUM_EVAL_INSTANCES eval instances...")
    t_gen = @elapsed begin
        eval_instances = [generate_er_instance(n_target, ER_P_EDGE;
            rng=MersenneTwister(SEED + n_target * 100 + i)) for i in 1:NUM_EVAL_INSTANCES]
    end
    eval_optimal = [maxcut_optimal(inst.costs) for inst in eval_instances]
    @printf("  Generated in %.1fs\n", t_gen)

    println("  Generating $NUM_HOMODIST_INSTANCES homodist instances...")
    t_gen2 = @elapsed begin
        homodist_pool = [generate_er_instance(n_target, ER_P_EDGE;
            rng=MersenneTwister(SEED + n_target * 200 + i)) for i in 1:NUM_HOMODIST_INSTANCES]
    end
    @printf("  Generated in %.1fs\n", t_gen2)

    for p in P_VALUES
        println("\n  p=$p:")

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
            med_γ = γs_t .* π
            med_β = βs_t .* π
        end

        println("    Evaluating Transfer...")
        t_tr = @elapsed begin
            transfer_ratios = map(enumerate(eval_instances)) do (i, inst)
                qaoa_expectation(inst.costs, n_target, med_γ, med_β) / eval_optimal[i]
            end
        end
        @printf("    Transfer:    mean=%.4f  std=%.4f  (%.1fs)\n",
            mean(transfer_ratios), std(transfer_ratios), t_tr)

        # === PAPERPROXY ===
        println("    Evaluating PaperProxy...")
        t_pp = @elapsed begin
            pp_ratios = map(enumerate(eval_instances)) do (i, inst)
                proxy = PaperProxy(inst.num_edges, n_target, ER_P_EDGE)
                homodist = cpu_compute_homodist(proxy)
                P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, n_target, p)
                qaoa_expectation(inst.costs, n_target, best_γ, best_β) / eval_optimal[i]
            end
        end
        @printf("    PaperProxy:  mean=%.4f  std=%.4f  (%.1fs)\n",
            mean(pp_ratios), std(pp_ratios), t_pp)

        # === SAMPLED HOMODIST ===
        println("    Computing sampled homodist (S=$SAMPLES_PER_COST)...")
        max_edges = maximum(inst.num_edges for inst in homodist_pool)
        t_sample = @elapsed begin
            sampled_homodists = map(homodist_pool) do inst
                get_homogeneous_distribution_from_costs_sampled(
                    inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
                    max_num_edges=max_edges,
                    rng=MersenneTwister(SEED + hash(inst.num_edges) + n_target))
            end
            sampled_avg = average_distributions(sampled_homodists)
        end
        P_emp = compute_empirical_P(homodist_pool, max_edges)
        m_hd = size(sampled_avg, 1) - 1
        P_vals = P_emp[1:min(m_hd+1, length(P_emp))]
        if length(P_vals) < m_hd + 1
            P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
        end
        @printf("    Sampled homodist computed in %.1fs\n", t_sample)

        println("    Optimizing via proxy...")
        t_opt = @elapsed begin
            best_γ_s, best_β_s, _ = optimize_via_proxy(sampled_avg, P_vals, n_target, p)
        end
        @printf("    Proxy optimization in %.1fs\n", t_opt)

        println("    Evaluating on target instances...")
        t_eval = @elapsed begin
            samp_ratios = map(enumerate(eval_instances)) do (i, inst)
                qaoa_expectation(inst.costs, n_target, best_γ_s, best_β_s) / eval_optimal[i]
            end
        end
        @printf("    SampN+EmpP:  mean=%.4f  std=%.4f  (total %.1fs)\n",
            mean(samp_ratios), std(samp_ratios), t_sample + t_opt + t_eval)
    end
end

println("\nDone!")
