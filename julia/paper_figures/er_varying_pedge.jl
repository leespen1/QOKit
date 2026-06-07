#=
er_varying_pedge.jl — ER with p_edge != 0.5: proxy and baseline comparison.

PaperProxy's analytical formula supports arbitrary p_edge. Sparse ER graphs
have much lower random baselines (0.51-0.60 vs 0.75 for ER(0.5)), so proxy
improvements are more impactful. This is a natural testbed for the proxy
heuristic.

Methods compared:
  (a) Random baseline (mean/max costs)
  (b) Transfer from small source graphs
  (c) PaperProxy with correct p_edge
  (d) SampN+EmpP (sampling-based homodist + empirical P)

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf
using Distributions: Binomial, pdf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const P_EDGES = [0.1, 0.2, 0.3, 0.5, 0.7]
const N_VALUES = [12, 14, 16, 18]
const N_SOURCE = 9
const NUM_HOMODIST_INSTANCES = 20
const NUM_EVAL_INSTANCES = 10
const SAMPLES_PER_COST = 20

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
println("ER with Varying p_edge: Proxy and Baseline Comparison")
println("=" ^ 80)

methods = ["Random", "Transfer", "PaperProxy", "SampN+EmpP"]

# Store all results: results[p_edge][method][n][p] = Vector{Float64}
all_results = Dict{Float64, Dict{String, Dict{Int, Dict{Int, Vector{Float64}}}}}()
for pe in P_EDGES
    all_results[pe] = Dict(m => Dict{Int, Dict{Int, Vector{Float64}}}() for m in methods)
end

for p_edge in P_EDGES
    println("\n" * "=" ^ 60)
    println("  p_edge = $p_edge")
    println("=" ^ 60)

    # Generate source instances for transfer
    src_rng = MersenneTwister(SEED + hash(p_edge) + 888)
    source_insts = [generate_er_instance(N_SOURCE, p_edge; rng=src_rng) for _ in 1:NUM_EVAL_INSTANCES]

    for n_target in N_VALUES
        # Skip very sparse at large n (too few edges for meaningful MaxCut)
        if p_edge <= 0.1 && n_target > 16
            continue
        end

        println("\n  --- n=$n_target, p_edge=$p_edge ---")

        # Generate eval and homodist instances
        eval_rng = MersenneTwister(SEED + hash(p_edge) + n_target * 100)
        eval_insts = [generate_er_instance(n_target, p_edge; rng=eval_rng) for _ in 1:NUM_EVAL_INSTANCES]
        eval_optimal = [maxcut_optimal(inst.costs) for inst in eval_insts]

        hd_rng = MersenneTwister(SEED + hash(p_edge) + n_target * 200)
        hd_insts = [generate_er_instance(n_target, p_edge; rng=hd_rng) for _ in 1:NUM_HOMODIST_INSTANCES]

        for p in P_VALUES
            println("  p=$p:")

            # === RANDOM BASELINE ===
            random_ratios = [mean(inst.costs) / eval_optimal[i] for (i, inst) in enumerate(eval_insts)]
            if !haskey(all_results[p_edge]["Random"], n_target)
                all_results[p_edge]["Random"][n_target] = Dict{Int, Vector{Float64}}()
            end
            all_results[p_edge]["Random"][n_target][p] = random_ratios
            print_stats("Random baseline", random_ratios)

            # === TRANSFER ===
            t_γ, t_β = transfer_optimize(source_insts, N_SOURCE, p)
            transfer_ratios = [qaoa_expectation(inst.costs, n_target, t_γ, t_β) / eval_optimal[i]
                               for (i, inst) in enumerate(eval_insts)]
            if !haskey(all_results[p_edge]["Transfer"], n_target)
                all_results[p_edge]["Transfer"][n_target] = Dict{Int, Vector{Float64}}()
            end
            all_results[p_edge]["Transfer"][n_target][p] = transfer_ratios
            print_stats("Transfer", transfer_ratios)

            # === PAPERPROXY (with correct p_edge) ===
            pp_ratios = map(enumerate(eval_insts)) do (i, inst)
                proxy = PaperProxy(inst.num_edges, n_target, p_edge)
                homodist = cpu_compute_homodist(proxy)
                P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, n_target, p)
                qaoa_expectation(inst.costs, n_target, best_γ, best_β) / eval_optimal[i]
            end
            if !haskey(all_results[p_edge]["PaperProxy"], n_target)
                all_results[p_edge]["PaperProxy"][n_target] = Dict{Int, Vector{Float64}}()
            end
            all_results[p_edge]["PaperProxy"][n_target][p] = pp_ratios
            print_stats("PaperProxy(p=$p_edge)", pp_ratios)

            # === SAMPLED HOMODIST ===
            max_edges = maximum(inst.num_edges for inst in hd_insts)
            t_sample = @elapsed begin
                sampled_homodists = map(hd_insts) do inst
                    get_homogeneous_distribution_from_costs_sampled(
                        inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
                        max_num_edges=max_edges,
                        rng=MersenneTwister(SEED + hash(p_edge) + hash(inst.num_edges) + n_target))
                end
                sampled_avg = average_distributions(sampled_homodists)
            end
            P_emp = compute_empirical_P(hd_insts, max_edges)
            m_hd = size(sampled_avg, 1) - 1
            P_vals = P_emp[1:min(m_hd+1, length(P_emp))]
            if length(P_vals) < m_hd + 1
                P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
            end

            best_γ_s, best_β_s, _ = optimize_via_proxy(sampled_avg, P_vals, n_target, p)
            samp_ratios = [qaoa_expectation(inst.costs, n_target, best_γ_s, best_β_s) / eval_optimal[i]
                           for (i, inst) in enumerate(eval_insts)]
            if !haskey(all_results[p_edge]["SampN+EmpP"], n_target)
                all_results[p_edge]["SampN+EmpP"][n_target] = Dict{Int, Vector{Float64}}()
            end
            all_results[p_edge]["SampN+EmpP"][n_target][p] = samp_ratios
            print_stats("SampN+EmpP ($(round(t_sample, digits=2))s)", samp_ratios)
        end
    end
end

#==============================================================================#
#   SUMMARY TABLES                                                             #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY TABLES")
println("=" ^ 80)

for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-6s %-4s", "p_edge", "n")
    for m in methods; @printf("  %-16s", m); end
    println()
    println("  " * "-" ^ (12 + 18 * length(methods)))
    for p_edge in P_EDGES
        for n in N_VALUES
            if p_edge <= 0.1 && n > 16; continue; end
            @printf("  %-6.1f %-4d", p_edge, n)
            for m in methods
                if haskey(all_results[p_edge][m], n) && haskey(all_results[p_edge][m][n], p)
                    vals = all_results[p_edge][m][n][p]
                    @printf("  %.4f±%.4f  ", mean(vals), std(vals))
                else
                    @printf("  %-16s", "—")
                end
            end
            println()
        end
    end
end

# Compute "improvement above random" for each method
println("\n\nIMPROVEMENT ABOVE RANDOM (method - random):")
for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-6s %-4s", "p_edge", "n")
    for m in methods[2:end]; @printf("  %-16s", m); end
    println()
    println("  " * "-" ^ (12 + 18 * (length(methods)-1)))
    for p_edge in P_EDGES
        for n in N_VALUES
            if p_edge <= 0.1 && n > 16; continue; end
            if !haskey(all_results[p_edge]["Random"], n); continue; end
            if !haskey(all_results[p_edge]["Random"][n], p); continue; end
            random_mean = mean(all_results[p_edge]["Random"][n][p])
            @printf("  %-6.1f %-4d", p_edge, n)
            for m in methods[2:end]
                if haskey(all_results[p_edge][m], n) && haskey(all_results[p_edge][m][n], p)
                    val_mean = mean(all_results[p_edge][m][n][p])
                    @printf("  %+.4f          ", val_mean - random_mean)
                else
                    @printf("  %-16s", "—")
                end
            end
            println()
        end
    end
end

println("\nDone!")
