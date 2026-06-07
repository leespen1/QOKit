#=
er_sparse_high_depth.jl — ER(0.2, 0.3) at p=1,3,5 and n=14-20.

On ER(0.5), proxy overfits at p=5. Does the same happen on sparse ER?
PaperProxy dominated at p=3 on sparse ER — if it persists at p=5, that's
a much stronger story than ER(0.5).

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const P_EDGES = [0.2, 0.3, 0.5]     # include 0.5 for comparison
const N_VALUES = [14, 16, 18, 20]
const N_SOURCE = 9
const NUM_HOMODIST_INSTANCES = 10
const NUM_EVAL_INSTANCES = 5
const SAMPLES_PER_COST = 20

const P_VALUES = [1, 3, 5]
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const SEED = 42

#==============================================================================#
#                    HELPERS (same as before)                                   #
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
        K = GRID_SIZE_P1^2
        γ_matrix = zeros(K, 1); β_matrix = zeros(K, 1)
        idx = 0
        for γ in range(0.02, 2.0, length=GRID_SIZE_P1), β in range(0.01, π/2 - 0.01, length=GRID_SIZE_P1)
            idx += 1; γ_matrix[idx, 1] = γ / π; β_matrix[idx, 1] = β / π
        end
        Qs = QAOA_proxy_multi(homodist, γ_matrix, β_matrix; pi_units=true)
        exps = vec(expectation(Qs[end], P_vals, n))
        best_idx = argmax(exps)
        return [γ_matrix[best_idx, 1] * π], [β_matrix[best_idx, 1] * π]
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
        return γ_matrix[best_idx, :] .* π, β_matrix[best_idx, :] .* π
    end
end

function transfer_params(source_insts, n_source, p)
    if p == 1
        sp = map(source_insts) do inst
            best_exp = -Inf; best_γ = 0.0; best_β = 0.0
            for γ in range(0.02, 2.0, length=80), β in range(0.01, π/2, length=80)
                e = qaoa_expectation(inst.costs, n_source, [γ], [β])
                if e > best_exp; best_exp = e; best_γ = γ; best_β = β; end
            end
            (γ=best_γ, β=best_β)
        end
        return [median([s.γ for s in sp])], [median([s.β for s in sp])]
    else
        gs = GRID_SIZE_RAMP
        sp = map(source_insts) do inst
            best_params = (0.0,0.0,0.0,0.0); best_exp = -Inf
            for γ₁ in range(0.02,0.40,length=gs), γ_f in range(0.10,0.70,length=gs),
                β₁ in range(0.05,0.45,length=gs), β_f in range(0.01,0.25,length=gs)
                γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, p)
                exp_val = qaoa_expectation(inst.costs, n_source, γs_pi .* π, βs_pi .* π)
                if exp_val > best_exp; best_exp = exp_val; best_params = (γ₁, γ_f, β₁, β_f); end
            end
            best_params
        end
        med = Tuple(median([s[i] for s in sp]) for i in 1:4)
        γs, βs = linear_ramp(med..., p)
        return γs .* π, βs .* π
    end
end

#==============================================================================#
#   MAIN                                                                       #
#==============================================================================#

println("=" ^ 80)
println("Sparse ER High-Depth Investigation: p=1,3,5")
println("=" ^ 80)

methods = ["Random", "Transfer", "PaperProxy", "SampN+EmpP"]

# results[p_edge][n][p][method] = mean_ratio
results = Dict{Float64, Dict{Int, Dict{Int, Dict{String, Float64}}}}()

for p_edge in P_EDGES
    results[p_edge] = Dict{Int, Dict{Int, Dict{String, Float64}}}()

    src_rng = MersenneTwister(SEED + hash(p_edge) + 888)
    source_insts = [generate_er_instance(N_SOURCE, p_edge; rng=src_rng) for _ in 1:5]

    for n in N_VALUES
        results[p_edge][n] = Dict{Int, Dict{String, Float64}}()

        eval_rng = MersenneTwister(SEED + hash(p_edge) + n * 100)
        eval_insts = [generate_er_instance(n, p_edge; rng=eval_rng) for _ in 1:NUM_EVAL_INSTANCES]
        eval_opt = [maxcut_optimal(inst.costs) for inst in eval_insts]

        hd_rng = MersenneTwister(SEED + hash(p_edge) + n * 200)
        hd_insts = [generate_er_instance(n, p_edge; rng=hd_rng) for _ in 1:NUM_HOMODIST_INSTANCES]

        max_edges = maximum(inst.num_edges for inst in hd_insts)
        sampled_homodists = map(hd_insts) do inst
            get_homogeneous_distribution_from_costs_sampled(
                inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
                max_num_edges=max_edges,
                rng=MersenneTwister(SEED + hash(p_edge) + hash(inst.num_edges) + n))
        end
        sampled_avg = average_distributions(sampled_homodists)
        P_emp = compute_empirical_P(hd_insts, max_edges)
        m_hd = size(sampled_avg, 1) - 1
        P_v = P_emp[1:min(m_hd+1, length(P_emp))]
        if length(P_v) < m_hd + 1; P_v = vcat(P_v, zeros(m_hd + 1 - length(P_v))); end

        for p in P_VALUES
            results[p_edge][n][p] = Dict{String, Float64}()

            # Random
            r = mean(mean(inst.costs) / eval_opt[i] for (i, inst) in enumerate(eval_insts))
            results[p_edge][n][p]["Random"] = r

            # Transfer
            tγ, tβ = transfer_params(source_insts, N_SOURCE, p)
            t = mean(qaoa_expectation(inst.costs, n, tγ, tβ) / eval_opt[i]
                     for (i, inst) in enumerate(eval_insts))
            results[p_edge][n][p]["Transfer"] = t

            # PaperProxy
            pp_vals = map(enumerate(eval_insts)) do (i, inst)
                proxy = PaperProxy(inst.num_edges, n, p_edge)
                hd = cpu_compute_homodist(proxy)
                Pv = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                bγ, bβ = optimize_via_proxy(hd, Pv, n, p)
                qaoa_expectation(inst.costs, n, bγ, bβ) / eval_opt[i]
            end
            results[p_edge][n][p]["PaperProxy"] = mean(pp_vals)

            # SampN+EmpP
            sγ, sβ = optimize_via_proxy(sampled_avg, P_v, n, p)
            s = mean(qaoa_expectation(inst.costs, n, sγ, sβ) / eval_opt[i]
                     for (i, inst) in enumerate(eval_insts))
            results[p_edge][n][p]["SampN+EmpP"] = s

            @printf("  p_e=%.1f n=%-2d p=%d  Rand=%.3f  Transf=%.3f  PP=%.3f  SampN=%.3f  best=%s\n",
                p_edge, n, p, r, t,
                results[p_edge][n][p]["PaperProxy"],
                s,
                argmax(Dict("PP" => results[p_edge][n][p]["PaperProxy"],
                           "SampN" => s, "Transfer" => t)))
        end
    end
end

# Summary table
println("\n" * "=" ^ 80)
println("SUMMARY: Best proxy method by (p_edge, p)")
println("=" ^ 80)
println()
@printf("%-6s %-4s %-4s  %-8s %-8s %-8s %-8s  %-12s %-12s\n",
    "p_edge", "n", "p", "Random", "Transfer", "PP", "SampN", "PP-Trans", "SampN-Trans")
println("-" ^ 90)
for p_edge in P_EDGES
    for n in N_VALUES
        for p in P_VALUES
            r = results[p_edge][n][p]
            @printf("%-6.1f %-4d %-4d  %.4f   %.4f   %.4f   %.4f   %+.4f      %+.4f\n",
                p_edge, n, p,
                r["Random"], r["Transfer"], r["PaperProxy"], r["SampN+EmpP"],
                r["PaperProxy"] - r["Transfer"],
                r["SampN+EmpP"] - r["Transfer"])
        end
    end
    println()
end

println("\nDone!")
