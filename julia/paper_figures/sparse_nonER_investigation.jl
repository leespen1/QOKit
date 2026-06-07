#=
sparse_nonER_investigation.jl — Test SampN+EmpP on sparse non-ER graphs.

Prior results showed SampN+EmpP degrades on non-ER at n=14-18, but those used
relatively dense BA(m=2) and WS(k=4). Sparse graphs have lower random baselines,
potentially giving more room for proxy improvement.

Test BA(m=1) and WS(k=2, p_rewire=0.3) which produce sparser graphs.
Also compare against dense versions for reference.

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VALUES = [12, 14, 16, 18]
const N_SOURCE = 9
const NUM_HOMODIST_INSTANCES = 20
const NUM_EVAL_INSTANCES = 10
const SAMPLES_PER_COST = 20

const P_VALUES = [1, 3]
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const SEED = 42

# Graph configs to test (name, generator_fn, params)
struct GraphConfig
    name::String
    type::String  # "BA" or "WS"
    # BA params
    m_attach::Int
    # WS params
    k::Int
    p_rewire::Float64
end

const GRAPH_CONFIGS = [
    GraphConfig("BA(m=1)", "BA", 1, 0, 0.0),
    GraphConfig("BA(m=2)", "BA", 2, 0, 0.0),
    GraphConfig("BA(m=3)", "BA", 3, 0, 0.0),
    GraphConfig("WS(k=2)", "WS", 0, 2, 0.3),
    GraphConfig("WS(k=4)", "WS", 0, 4, 0.3),
    GraphConfig("WS(k=6)", "WS", 0, 6, 0.3),
]

#==============================================================================#
#                    HELPERS                                                    #
#==============================================================================#

function generate_config_instance(cfg::GraphConfig, n::Int; rng=Random.default_rng())
    if cfg.type == "BA"
        return generate_ba_instance(n, cfg.m_attach; rng)
    else
        return generate_ws_instance(n, cfg.k, cfg.p_rewire; rng)
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

function transfer_optimize(source_instances, n_source, p)
    if p == 1
        src_params = map(source_instances) do inst
            best_exp = -Inf; best_γ = 0.0; best_β = 0.0
            for γ in range(0.02, 2.0, length=80), β in range(0.01, π/2, length=80)
                e = qaoa_expectation(inst.costs, n_source, [γ], [β])
                if e > best_exp; best_exp = e; best_γ = γ; best_β = β; end
            end
            (γ=best_γ, β=best_β)
        end
        return [median([sp.γ for sp in src_params])], [median([sp.β for sp in src_params])]
    else
        gs = GRID_SIZE_RAMP
        src_ramp = map(source_instances) do inst
            best_params = (0.0, 0.0, 0.0, 0.0); best_exp = -Inf
            for γ₁ in range(0.02, 0.40, length=gs),
                γ_f in range(0.10, 0.70, length=gs),
                β₁ in range(0.05, 0.45, length=gs),
                β_f in range(0.01, 0.25, length=gs)
                γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, p)
                exp_val = qaoa_expectation(inst.costs, n_source, γs_pi .* π, βs_pi .* π)
                if exp_val > best_exp; best_exp = exp_val; best_params = (γ₁, γ_f, β₁, β_f); end
            end
            best_params
        end
        med_ramp = Tuple(median([sp[i] for sp in src_ramp]) for i in 1:4)
        γs, βs = linear_ramp(med_ramp..., p)
        return γs .* π, βs .* π
    end
end

function print_stats(label, values)
    @printf("    %-50s  mean=%.4f  std=%.4f\n", label, mean(values), std(values))
end

#==============================================================================#
#   MAIN EXPERIMENT                                                            #
#==============================================================================#

println("=" ^ 80)
println("Sparse Non-ER Investigation: SampN+EmpP on BA/WS with varying density")
println("=" ^ 80)

methods = ["Random", "Transfer", "PaperProxy", "SampN+EmpP"]

# results[config_name][method][n][p] = Vector{Float64}
all_results = Dict{String, Dict{String, Dict{Int, Dict{Int, Vector{Float64}}}}}()
for cfg in GRAPH_CONFIGS
    all_results[cfg.name] = Dict(m => Dict{Int, Dict{Int, Vector{Float64}}}() for m in methods)
end

for cfg in GRAPH_CONFIGS
    # BA(m=1) doesn't work at n_source=9 with m_attach=1 (too small)
    n_source = max(N_SOURCE, cfg.m_attach + 2)

    for n_target in N_VALUES
        # Skip n<k+1 for WS
        if cfg.type == "WS" && n_target <= cfg.k; continue; end

        println("\n" * "=" ^ 40)
        println("  $(cfg.name) (n=$n_target)")
        println("=" ^ 40)

        # Generate instances
        eval_rng = MersenneTwister(SEED + hash(cfg.name) + n_target * 100)
        eval_insts = [generate_config_instance(cfg, n_target; rng=eval_rng) for _ in 1:NUM_EVAL_INSTANCES]
        eval_optimal = [maxcut_optimal(inst.costs) for inst in eval_insts]

        hd_rng = MersenneTwister(SEED + hash(cfg.name) + n_target * 200)
        hd_insts = [generate_config_instance(cfg, n_target; rng=hd_rng) for _ in 1:NUM_HOMODIST_INSTANCES]

        src_rng = MersenneTwister(SEED + hash(cfg.name) + 888)
        source_insts = [generate_config_instance(cfg, n_source; rng=src_rng) for _ in 1:NUM_EVAL_INSTANCES]

        # Stats
        m_avg = mean(inst.num_edges for inst in eval_insts)
        p_eff_avg = mean(effective_edge_probability(inst.num_edges, n_target) for inst in eval_insts)
        @printf("  avg edges=%.1f, p_eff=%.3f\n", m_avg, p_eff_avg)

        # Compute sampled homodist
        max_edges = maximum(inst.num_edges for inst in hd_insts)
        t_hd = @elapsed begin
            sampled_homodists = map(hd_insts) do inst
                get_homogeneous_distribution_from_costs_sampled(
                    inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
                    max_num_edges=max_edges,
                    rng=MersenneTwister(SEED + hash(cfg.name) + hash(inst.num_edges) + n_target))
            end
            sampled_avg = average_distributions(sampled_homodists)
        end
        P_emp = compute_empirical_P(hd_insts, max_edges)
        m_hd = size(sampled_avg, 1) - 1
        P_vals = P_emp[1:min(m_hd+1, length(P_emp))]
        if length(P_vals) < m_hd + 1
            P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
        end

        for p in P_VALUES
            println("  p=$p:")

            # === RANDOM ===
            random_ratios = [mean(inst.costs) / eval_optimal[i] for (i, inst) in enumerate(eval_insts)]
            get!(all_results[cfg.name]["Random"], n_target, Dict{Int,Vector{Float64}}())[p] = random_ratios
            print_stats("Random baseline", random_ratios)

            # === TRANSFER ===
            t_γ, t_β = transfer_optimize(source_insts, n_source, p)
            transfer_ratios = [qaoa_expectation(inst.costs, n_target, t_γ, t_β) / eval_optimal[i]
                               for (i, inst) in enumerate(eval_insts)]
            get!(all_results[cfg.name]["Transfer"], n_target, Dict{Int,Vector{Float64}}())[p] = transfer_ratios
            print_stats("Transfer", transfer_ratios)

            # === PAPERPROXY ===
            pp_ratios = map(enumerate(eval_insts)) do (i, inst)
                p_eff = effective_edge_probability(inst.num_edges, n_target)
                proxy = PaperProxy(inst.num_edges, n_target, p_eff)
                hd = cpu_compute_homodist(proxy)
                Pv = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                best_γ, best_β, _ = optimize_via_proxy(hd, Pv, n_target, p)
                qaoa_expectation(inst.costs, n_target, best_γ, best_β) / eval_optimal[i]
            end
            get!(all_results[cfg.name]["PaperProxy"], n_target, Dict{Int,Vector{Float64}}())[p] = pp_ratios
            print_stats("PaperProxy(p_eff)", pp_ratios)

            # === SAMPN+EMPP ===
            best_γ_s, best_β_s, _ = optimize_via_proxy(sampled_avg, P_vals, n_target, p)
            samp_ratios = [qaoa_expectation(inst.costs, n_target, best_γ_s, best_β_s) / eval_optimal[i]
                           for (i, inst) in enumerate(eval_insts)]
            get!(all_results[cfg.name]["SampN+EmpP"], n_target, Dict{Int,Vector{Float64}}())[p] = samp_ratios
            print_stats("SampN+EmpP ($(round(t_hd, digits=2))s)", samp_ratios)
        end
    end
end

#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-10s %-4s %-8s", "Config", "n", "p_eff")
    for m in methods; @printf("  %-16s", m); end
    println()
    println("  " * "-" ^ (24 + 18 * length(methods)))
    for cfg in GRAPH_CONFIGS
        for n in N_VALUES
            if cfg.type == "WS" && n <= cfg.k; continue; end
            if !haskey(all_results[cfg.name]["Random"], n); continue; end
            if !haskey(all_results[cfg.name]["Random"][n], p); continue; end

            # Compute average p_eff from results
            eval_rng = MersenneTwister(SEED + hash(cfg.name) + n * 100)
            insts = [generate_config_instance(cfg, n; rng=eval_rng) for _ in 1:3]
            p_eff = mean(effective_edge_probability(inst.num_edges, n) for inst in insts)

            @printf("  %-10s %-4d %-8.3f", cfg.name, n, p_eff)
            for m in methods
                if haskey(all_results[cfg.name][m], n) && haskey(all_results[cfg.name][m][n], p)
                    vals = all_results[cfg.name][m][n][p]
                    @printf("  %.4f±%.4f  ", mean(vals), std(vals))
                else
                    @printf("  %-16s", "—")
                end
            end
            println()
        end
    end
end

# Compute "gap above random" for SampN+EmpP vs Transfer
println("\n\nSampN+EmpP - Transfer DIFFERENCE (positive = SampN beats Transfer):")
for p in P_VALUES
    println("\n  p=$p:")
    for cfg in GRAPH_CONFIGS
        for n in N_VALUES
            if !haskey(all_results[cfg.name]["SampN+EmpP"], n); continue; end
            if !haskey(all_results[cfg.name]["SampN+EmpP"][n], p); continue; end
            s = mean(all_results[cfg.name]["SampN+EmpP"][n][p])
            t = mean(all_results[cfg.name]["Transfer"][n][p])
            @printf("  %-10s n=%-4d  SampN=%.4f  Transfer=%.4f  diff=%+.4f\n",
                cfg.name, n, s, t, s - t)
        end
    end
end

println("\nDone!")
