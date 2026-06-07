#=
high_depth_large_n_er.jl — Test SampN+EmpP at higher depths (p=5) on large ER graphs.

The headline result (SampN+EmpP beats Transfer by 3.6% at n=22,p=3) only covers
p=1 and p=3. This experiment extends to p=5 to verify the advantage persists
at higher depth with linear ramp schedules.

At large n, ramp optimization uses a 4D grid (γ₁,γ_f,β₁,β_f). For p=5, the
proxy evaluation is fast but QAOA eval is O(2^n * p), so we keep n<=20 for p=5.

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VALUES = [14, 16, 18, 20]
const N_SOURCE = 9
const NUM_HOMODIST_INSTANCES = 10
const NUM_EVAL_INSTANCES = 5
const SAMPLES_PER_COST = 20

const ER_P_EDGE = 0.5

const P_VALUES = [1, 3, 5]
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
println("SampN+EmpP at Higher Depths on Large ER Graphs")
println("n = $N_VALUES, p = $P_VALUES")
println("=" ^ 80)

# Source instances for transfer
source_instances = [generate_er_instance(N_SOURCE, ER_P_EDGE;
    rng=MersenneTwister(SEED + 888 + i)) for i in 1:10]

# Store results for summary
all_results = Dict{Tuple{Int,Int}, Dict{String, Vector{Float64}}}()

for n_target in N_VALUES
    println("\n" * "=" ^ 40)
    @printf("  n = %d  (2^n = %d)\n", n_target, 1 << n_target)
    println("=" ^ 40)

    # Generate instances
    println("  Generating instances...")
    eval_instances = [generate_er_instance(n_target, ER_P_EDGE;
        rng=MersenneTwister(SEED + n_target * 100 + i)) for i in 1:NUM_EVAL_INSTANCES]
    eval_optimal = [maxcut_optimal(inst.costs) for inst in eval_instances]

    homodist_pool = [generate_er_instance(n_target, ER_P_EDGE;
        rng=MersenneTwister(SEED + n_target * 200 + i)) for i in 1:NUM_HOMODIST_INSTANCES]

    # Precompute sampled homodist (shared across depths)
    max_edges = maximum(inst.num_edges for inst in homodist_pool)
    println("  Computing sampled homodist (S=$SAMPLES_PER_COST)...")
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
    P_vals_samp = P_emp[1:min(m_hd+1, length(P_emp))]
    if length(P_vals_samp) < m_hd + 1
        P_vals_samp = vcat(P_vals_samp, zeros(m_hd + 1 - length(P_vals_samp)))
    end
    @printf("  Sampled homodist in %.1fs\n", t_sample)

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

        t_tr = @elapsed begin
            transfer_ratios = map(enumerate(eval_instances)) do (i, inst)
                qaoa_expectation(inst.costs, n_target, med_γ, med_β) / eval_optimal[i]
            end
        end
        @printf("    Transfer:    mean=%.4f  std=%.4f  (%.1fs)\n",
            mean(transfer_ratios), std(transfer_ratios), t_tr)

        # === PAPERPROXY ===
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
        t_opt = @elapsed begin
            best_γ_s, best_β_s, _ = optimize_via_proxy(sampled_avg, P_vals_samp, n_target, p)
        end
        t_eval = @elapsed begin
            samp_ratios = map(enumerate(eval_instances)) do (i, inst)
                qaoa_expectation(inst.costs, n_target, best_γ_s, best_β_s) / eval_optimal[i]
            end
        end
        @printf("    SampN+EmpP:  mean=%.4f  std=%.4f  (opt %.1fs, eval %.1fs)\n",
            mean(samp_ratios), std(samp_ratios), t_opt, t_eval)

        all_results[(n_target, p)] = Dict(
            "Transfer" => transfer_ratios,
            "PaperProxy" => pp_ratios,
            "SampN+EmpP" => samp_ratios
        )
    end
end

#==============================================================================#
#   SUMMARY TABLE                                                              #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY — Mean Approximation Ratios (ER, p=0.5)")
println("=" ^ 80)

@printf("%-6s %-4s %-12s %-12s %-12s %-10s\n", "n", "p", "Transfer", "PaperProxy", "SampN+EmpP", "Gap(S-T)")
println("-" ^ 60)
for n in N_VALUES
    for p in P_VALUES
        key = (n, p)
        haskey(all_results, key) || continue
        r = all_results[key]
        t_mean = mean(r["Transfer"])
        pp_mean = mean(r["PaperProxy"])
        s_mean = mean(r["SampN+EmpP"])
        gap = s_mean - t_mean
        marker = gap > 0.005 ? " ***" : gap > 0 ? " *" : ""
        @printf("%-6d %-4d %-12.4f %-12.4f %-12.4f %+.4f%s\n",
            n, p, t_mean, pp_mean, s_mean, gap, marker)
    end
end

println("\nDone!")
