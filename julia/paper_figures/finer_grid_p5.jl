#=
finer_grid_p5.jl — Test whether a finer ramp grid helps SampN+EmpP at p=5.

The high-depth experiment showed Transfer beating SampN+EmpP at p=5. One
hypothesis is that the 10^4 ramp grid is too coarse. Test 15^4 = 50625 and
20^4 = 160000 grid sizes.

Quick experiment: only n=18 on ER, p=5.

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_TARGET = 18
const N_SOURCE = 9
const NUM_HOMODIST_INSTANCES = 10
const NUM_EVAL_INSTANCES = 5
const SAMPLES_PER_COST = 20
const ER_P_EDGE = 0.5
const P_DEPTH = 5
const GRID_SIZES = [10, 15, 20]  # per dimension → 10^4, 15^4, 20^4 total
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

function optimize_via_proxy_gs(homodist, P_vals, n, p, gs)
    K = gs^4
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

#==============================================================================#
#   MAIN EXPERIMENT                                                            #
#==============================================================================#

println("=" ^ 80)
println("Finer Grid for SampN+EmpP at p=$P_DEPTH, n=$N_TARGET (ER)")
println("=" ^ 80)

# Instances
eval_instances = [generate_er_instance(N_TARGET, ER_P_EDGE;
    rng=MersenneTwister(SEED + N_TARGET * 100 + i)) for i in 1:NUM_EVAL_INSTANCES]
eval_optimal = [maxcut_optimal(inst.costs) for inst in eval_instances]

homodist_pool = [generate_er_instance(N_TARGET, ER_P_EDGE;
    rng=MersenneTwister(SEED + N_TARGET * 200 + i)) for i in 1:NUM_HOMODIST_INSTANCES]

# Sampled homodist
max_edges = maximum(inst.num_edges for inst in homodist_pool)
sampled_homodists = map(homodist_pool) do inst
    get_homogeneous_distribution_from_costs_sampled(
        inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
        max_num_edges=max_edges,
        rng=MersenneTwister(SEED + hash(inst.num_edges) + N_TARGET))
end
sampled_avg = average_distributions(sampled_homodists)
P_emp = compute_empirical_P(homodist_pool, max_edges)
m_hd = size(sampled_avg, 1) - 1
P_vals = P_emp[1:min(m_hd+1, length(P_emp))]
if length(P_vals) < m_hd + 1
    P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
end

# Transfer baseline
source_instances = [generate_er_instance(N_SOURCE, ER_P_EDGE;
    rng=MersenneTwister(SEED + 888 + i)) for i in 1:10]
src_ramp = map(source_instances) do inst
    best_params = (0.0, 0.0, 0.0, 0.0); best_exp = -Inf; gs = 10
    for γ₁ in range(0.02, 0.40, length=gs),
        γ_f in range(0.10, 0.70, length=gs),
        β₁ in range(0.05, 0.45, length=gs),
        β_f in range(0.01, 0.25, length=gs)
        γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, P_DEPTH)
        exp_val = qaoa_expectation(inst.costs, N_SOURCE, γs_pi .* π, βs_pi .* π)
        if exp_val > best_exp; best_exp = exp_val; best_params = (γ₁, γ_f, β₁, β_f); end
    end
    best_params
end
med_ramp = Tuple(median([sp[i] for sp in src_ramp]) for i in 1:4)
γs_t, βs_t = linear_ramp(med_ramp..., P_DEPTH)

transfer_ratios = map(enumerate(eval_instances)) do (i, inst)
    qaoa_expectation(inst.costs, N_TARGET, γs_t .* π, βs_t .* π) / eval_optimal[i]
end
@printf("Transfer (gs=10):  mean=%.4f  std=%.4f\n", mean(transfer_ratios), std(transfer_ratios))

# Also test Transfer with finer grids on source
for gs_src in [15, 20]
    src_ramp_fine = map(source_instances) do inst
        best_params = (0.0, 0.0, 0.0, 0.0); best_exp = -Inf
        for γ₁ in range(0.02, 0.40, length=gs_src),
            γ_f in range(0.10, 0.70, length=gs_src),
            β₁ in range(0.05, 0.45, length=gs_src),
            β_f in range(0.01, 0.25, length=gs_src)
            γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, P_DEPTH)
            exp_val = qaoa_expectation(inst.costs, N_SOURCE, γs_pi .* π, βs_pi .* π)
            if exp_val > best_exp; best_exp = exp_val; best_params = (γ₁, γ_f, β₁, β_f); end
        end
        best_params
    end
    med_fine = Tuple(median([sp[i] for sp in src_ramp_fine]) for i in 1:4)
    γs_tf, βs_tf = linear_ramp(med_fine..., P_DEPTH)
    tr_fine = map(enumerate(eval_instances)) do (i, inst)
        qaoa_expectation(inst.costs, N_TARGET, γs_tf .* π, βs_tf .* π) / eval_optimal[i]
    end
    @printf("Transfer (gs=%d):  mean=%.4f  std=%.4f\n", gs_src, mean(tr_fine), std(tr_fine))
end

# SampN+EmpP with varying grid sizes
for gs in GRID_SIZES
    K = gs^4
    @printf("\nSampN+EmpP (gs=%d, K=%d):\n", gs, K)
    t = @elapsed begin
        best_γ, best_β, proxy_exp = optimize_via_proxy_gs(sampled_avg, P_vals, N_TARGET, P_DEPTH, gs)
    end
    ratios = map(enumerate(eval_instances)) do (i, inst)
        qaoa_expectation(inst.costs, N_TARGET, best_γ, best_β) / eval_optimal[i]
    end
    @printf("  mean=%.4f  std=%.4f  proxy_exp=%.4f  (%.1fs)\n",
        mean(ratios), std(ratios), proxy_exp, t)
end

println("\nDone!")
