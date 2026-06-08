#=
smallworld_investigation.jl — Investigate QAOA proxy performance on
Barabási-Albert and Watts-Strogatz graphs.

Questions addressed:
  1. Do P(c') and N(c';d,c) differ from Erdős-Rényi?
  2. Does parameter transfer (small → large graphs) work?
  3. Does the PaperProxy (with effective edge probability) work?
  4. How does PaperProxy compare to parameter transfer?

For fair comparison at p>1, both Transfer and Proxy use linear ramp schedules
optimized by the same grid search. At p=1, Transfer uses coordinate descent
and Proxy uses 2D grid search (both are thorough for 2 parameters).

Assumptions:
  - BA: m_attach=2, WS: k=4 p_rewire=0.3
  - PaperProxy on non-ER: p_eff = m / (n choose 2)
  - For ER: PaperProxy uses known p_edge=0.5

Started: 2026-03-17
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_SMALL = 9
const N_LARGE = 12
const NUM_INSTANCES = 20

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

const P_VALUES = [1, 2, 3]
const N_RESTARTS = 15            # For p=1 coordinate descent
const GRID_SIZE_P1 = 50          # p=1 grid size (50^2 = 2500 points)
const GRID_SIZE_RAMP = 10        # Linear ramp grid (10^4 = 10000 points)
const SEED = 42
const N_DIST = 10                # n for distribution comparison


#==============================================================================#
#                    HELPER FUNCTIONS                                           #
#==============================================================================#

"""Optimize QAOA parameters by coordinate descent (for p=1 mainly)."""
function optimize_qaoa_cd(costs, n, p, n_restarts; rng=Random.default_rng())
    best_exp = -Inf
    best_γs = zeros(p)
    best_βs = zeros(p)

    for _ in 1:n_restarts
        γs = rand(rng, p) .* 1.6
        βs = rand(rng, p) .* (π/2)
        current_γs, current_βs = copy(γs), copy(βs)
        current_exp = qaoa_expectation(costs, n, current_γs, current_βs)

        for step_scale in [0.3, 0.15, 0.07, 0.03, 0.01]
            for param_idx in 1:(2p)
                for delta in [-2, -1, -0.5, 0.5, 1, 2] .* step_scale
                    trial_γs, trial_βs = copy(current_γs), copy(current_βs)
                    if param_idx <= p
                        trial_γs[param_idx] = max(0, trial_γs[param_idx] + delta)
                    else
                        trial_βs[param_idx-p] = clamp(trial_βs[param_idx-p] + delta, 0, π/2)
                    end
                    trial_exp = qaoa_expectation(costs, n, trial_γs, trial_βs)
                    if trial_exp > current_exp
                        current_exp = trial_exp
                        current_γs, current_βs = trial_γs, trial_βs
                    end
                end
            end
        end

        if current_exp > best_exp
            best_exp = current_exp
            best_γs, best_βs = copy(current_γs), copy(current_βs)
        end
    end
    return best_γs, best_βs, best_exp
end

"""Optimize QAOA linear ramp on REAL QAOA (not proxy) via grid search."""
function optimize_qaoa_linramp(costs, n, p; grid_size=GRID_SIZE_RAMP)
    γ₁_range = range(0.02, 0.40, length=grid_size)
    γ_f_range = range(0.10, 0.70, length=grid_size)
    β₁_range = range(0.05, 0.45, length=grid_size)
    β_f_range = range(0.01, 0.25, length=grid_size)

    best_exp = -Inf
    best_params = (0.0, 0.0, 0.0, 0.0)

    for γ₁ in γ₁_range, γ_f in γ_f_range, β₁ in β₁_range, β_f in β_f_range
        γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, p)
        γs_rad = γs_pi .* π
        βs_rad = βs_pi .* π
        exp_val = qaoa_expectation(costs, n, γs_rad, βs_rad)
        if exp_val > best_exp
            best_exp = exp_val
            best_params = (γ₁, γ_f, β₁, β_f)
        end
    end

    γs_pi, βs_pi = linear_ramp(best_params..., p)
    return best_params, γs_pi .* π, βs_pi .* π, best_exp
end

"""Optimize via proxy grid search."""
function optimize_via_proxy(homodist, P_vals, n, p)
    m = size(homodist, 1) - 1
    if p == 1
        γ_range = range(0.02, 2.0, length=GRID_SIZE_P1)
        β_range = range(0.01, π/2 - 0.01, length=GRID_SIZE_P1)
        K = GRID_SIZE_P1^2
        γ_matrix = zeros(K, 1)
        β_matrix = zeros(K, 1)
        idx = 0
        for γ in γ_range, β in β_range
            idx += 1
            γ_matrix[idx, 1] = γ / π
            β_matrix[idx, 1] = β / π
        end
        Qs = QAOA_proxy_multi(homodist, γ_matrix, β_matrix; pi_units=true)
        exps = vec(expectation(Qs[end], P_vals, n))
        best_idx = argmax(exps)
        return [γ_matrix[best_idx, 1] * π], [β_matrix[best_idx, 1] * π], exps[best_idx]
    else
        # Linear ramp with same grid as transfer_linramp
        gs = GRID_SIZE_RAMP
        γ₁_range = range(0.02, 0.40, length=gs)
        γ_f_range = range(0.10, 0.70, length=gs)
        β₁_range = range(0.05, 0.45, length=gs)
        β_f_range = range(0.01, 0.25, length=gs)

        K = gs^4
        γ_matrix = zeros(K, p)
        β_matrix = zeros(K, p)
        idx = 0
        for γ₁ in γ₁_range, γ_f in γ_f_range, β₁ in β₁_range, β_f in β_f_range
            idx += 1
            γs, βs = linear_ramp(γ₁, γ_f, β₁, β_f, p)
            γ_matrix[idx, :] .= γs
            β_matrix[idx, :] .= βs
        end

        Qs = QAOA_proxy_multi(homodist, γ_matrix, β_matrix; pi_units=true)
        exps = vec(expectation(Qs[end], P_vals, n))
        best_idx = argmax(exps)
        return γ_matrix[best_idx, :] .* π, β_matrix[best_idx, :] .* π, exps[best_idx]
    end
end

"""Generate instances for a given graph type."""
function generate_instances(graph_type::String, n::Int, num::Int; rng=Random.default_rng())
    map(1:num) do _
        if graph_type == "ER"
            generate_er_instance(n, ER_P_EDGE; rng)
        elseif graph_type == "BA"
            generate_ba_instance(n, BA_M_ATTACH; rng)
        elseif graph_type == "WS"
            generate_ws_instance(n, WS_K, WS_P_REWIRE; rng)
        else
            error("Unknown: $graph_type")
        end
    end
end

"""Compute PaperProxy for a graph instance."""
function make_proxy(graph_type::String, num_edges::Int, num_vertices::Int)
    p_e = graph_type == "ER" ? ER_P_EDGE : effective_edge_probability(num_edges, num_vertices)
    PaperProxy(num_edges, num_vertices, p_e)
end

function print_stats(label, values)
    @printf("  %-40s  mean=%.4f  std=%.4f  [%.4f, %.4f]\n",
        label, mean(values), std(values), minimum(values), maximum(values))
end


#==============================================================================#
#   QUESTION 1: Do P(c') and N(c';d,c) differ from ER?                        #
#==============================================================================#

println("=" ^ 80)
println("QUESTION 1: Distribution Comparison")
println("=" ^ 80)
println("n=$N_DIST, $NUM_INSTANCES instances per graph type\n")

graph_types = ["ER", "BA", "WS"]
dist_instances = Dict(gt => generate_instances(gt, N_DIST, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt))) for gt in graph_types)

println("--- Edge Count Statistics ---")
for gt in graph_types
    edges = [inst.num_edges for inst in dist_instances[gt]]
    print_stats("$gt edges", Float64.(edges))
end

# Empirical P(c')
max_edges_all = maximum(inst.num_edges for gt in graph_types for inst in dist_instances[gt])

println("\n--- Empirical P(c') Statistics ---")
empirical_P = Dict{String, Vector{Float64}}()
for gt in graph_types
    P_avg = zeros(max_edges_all + 1)
    for inst in dist_instances[gt]
        num_bs = 1 << N_DIST
        for x in 1:num_bs
            c = Int(inst.costs[x])
            if c <= max_edges_all
                P_avg[c + 1] += 1.0 / (num_bs * NUM_INSTANCES)
            end
        end
    end
    empirical_P[gt] = P_avg
    mean_cost = sum(c * P_avg[c+1] for c in 0:max_edges_all)
    var_cost = sum((c - mean_cost)^2 * P_avg[c+1] for c in 0:max_edges_all)
    @printf("  %-5s  mean_cost=%.2f  std_cost=%.2f\n", gt, mean_cost, sqrt(var_cost))
end

# Empirical N(c';d,c)
println("\n--- Computing Empirical N(c';d,c) ---")
empirical_homodists = Dict{String, Array{Float64, 3}}()
for gt in graph_types
    println("  $gt...")
    homodists = map(dist_instances[gt]) do inst
        get_homogeneous_distribution_from_costs_direct(
            inst.costs, inst.num_edges, inst.num_vertices;
            max_num_edges=max_edges_all)
    end
    empirical_homodists[gt] = average_distributions(homodists)
end

# PaperProxy
println("\n--- PaperProxy(p_eff) vs Empirical ---")
proxy_homodists = Dict{String, Array{Float64, 3}}()
for gt in graph_types
    avg_edges = round(Int, mean(inst.num_edges for inst in dist_instances[gt]))
    proxy = make_proxy(gt, avg_edges, N_DIST)
    p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(avg_edges, N_DIST)
    @printf("  %-5s  avg_edges=%d  p_eff=%.4f\n", gt, avg_edges, p_eff)
    proxy_homodists[gt] = cpu_compute_homodist(proxy)
end

# Pearson correlations
println("\n--- Pearson Correlations ---")
for gt in graph_types
    ph, eh = pad_to_match(proxy_homodists[gt], empirical_homodists[gt])
    corrs = get_pearson_correlation_coefficients(ph, eh)
    P_vals = empirical_P[gt]
    num_costs = min(length(corrs), length(P_vals))
    valid = [i for i in 1:num_costs if !isnan(corrs[i]) && P_vals[i] > 0.01]
    if !isempty(valid)
        weighted_corr = sum(corrs[i] * P_vals[i] for i in valid) / sum(P_vals[i] for i in valid)
        min_corr = minimum(corrs[i] for i in valid)
        @printf("  %-5s  P-weighted_corr=%.4f  min_corr(dominant)=%.4f  n_dominant=%d\n",
            gt, weighted_corr, min_corr, length(valid))
    end
end

# Also compute normalized MSE between proxy and empirical (not just Pearson)
println("\n--- Normalized MSE: Proxy vs Empirical ---")
for gt in graph_types
    ph, eh = pad_to_match(proxy_homodists[gt], empirical_homodists[gt])
    # Normalize each so that sum = 1
    ph_n = ph ./ max(sum(ph), 1e-10)
    eh_n = eh ./ max(sum(eh), 1e-10)
    mse = sum((ph_n .- eh_n).^2) / length(ph_n)
    @printf("  %-5s  normalized_MSE=%.6e\n", gt, mse)
end

# Homogeneity check
println("\n--- Homogeneity Check (CV of N(c';d,c) across instances) ---")
for gt in graph_types
    homodists = map(dist_instances[gt]) do inst
        get_homogeneous_distribution_from_costs_direct(
            inst.costs, inst.num_edges, inst.num_vertices;
            max_num_edges=max_edges_all)
    end
    m_arr, s_arr = distributions_mean_and_stddev(homodists)
    cv_vals = [s_arr[idx] / m_arr[idx] for idx in CartesianIndices(m_arr) if m_arr[idx] > 1.0]
    if !isempty(cv_vals)
        @printf("  %-5s  median_CV=%.4f  mean_CV=%.4f  (%d dominant entries)\n",
            gt, median(cv_vals), mean(cv_vals), length(cv_vals))
    end
end

# Cross-type MSE
println("\n--- Cross-Type Distribution MSE ---")
for gt1 in graph_types, gt2 in graph_types
    if gt1 >= gt2; continue; end
    h1, h2 = pad_to_match(empirical_homodists[gt1], empirical_homodists[gt2])
    h1n = h1 ./ max(sum(h1), 1e-10)
    h2n = h2 ./ max(sum(h2), 1e-10)
    @printf("  %-5s vs %-5s  normalized_MSE=%.6e\n", gt1, gt2,
        sum((h1n .- h2n).^2) / length(h1n))
end


#==============================================================================#
#   Q1 PLOTS                                                                   #
#==============================================================================#

println("\n--- Q1 Plots ---")

fig1 = Figure(size=(1000, 400))
gt_colors = Dict("ER" => :steelblue, "BA" => :coral, "WS" => :mediumseagreen)

ax1 = Axis(fig1[1, 1], xlabel="Cost c'", ylabel="P(c')",
    title="Empirical P(c')\n(n=$N_DIST, $NUM_INSTANCES instances)")
for gt in graph_types
    P = empirical_P[gt]
    lines!(ax1, 0:length(P)-1, P, label=gt, color=gt_colors[gt], linewidth=2)
    scatter!(ax1, 0:length(P)-1, P, color=gt_colors[gt], markersize=4)
end
axislegend(ax1, position=:rt)

ax2 = Axis(fig1[1, 2], xlabel="Cost c'", ylabel="Pearson Correlation",
    title="Correlation: Empirical vs PaperProxy(p_eff)")
for gt in graph_types
    ph, eh = pad_to_match(proxy_homodists[gt], empirical_homodists[gt])
    corrs = replace(get_pearson_correlation_coefficients(ph, eh), NaN => 0.0)
    scatterlines!(ax2, 0:length(corrs)-1, corrs, label=gt, color=gt_colors[gt], markersize=6)
end
hlines!(ax2, [0.9], color=:gray, linestyle=:dash, label="0.9")
axislegend(ax2, position=:lb)

save_figure(fig1, "smallworld_q1_distributions.png")

# Heatmaps
c_prime_show = max_edges_all ÷ 2
fig1b = Figure(size=(1200, 350 * length(graph_types)))
for (gi, gt) in enumerate(graph_types)
    emp = empirical_homodists[gt]
    ci = min(c_prime_show + 1, size(emp, 1))
    emp_s = emp[ci, :, :]
    nd, nc = size(emp_s)
    vmax = max(maximum(emp_s), 1e-10)

    ax_e = Axis(fig1b[gi, 1], xlabel="Cost c", ylabel="Dist d", title="$gt Empirical, c'=$(ci-1)")
    heatmap!(ax_e, 0:nc-1, 0:nd-1, emp_s, colorrange=(0, vmax), colormap=:viridis)

    ph, _ = pad_to_match(proxy_homodists[gt], emp)
    ps = ph[ci, 1:nd, 1:nc]
    ax_p = Axis(fig1b[gi, 2], xlabel="Cost c", ylabel="Dist d", title="$gt PaperProxy(p_eff), c'=$(ci-1)")
    hm = heatmap!(ax_p, 0:nc-1, 0:nd-1, ps, colorrange=(0, vmax), colormap=:viridis)
    Colorbar(fig1b[gi, 3], hm, label="N(c';d,c)")

    ds = emp_s .- ps
    dmax = max(maximum(abs.(ds)), 1e-10)
    ax_d = Axis(fig1b[gi, 4], xlabel="Cost c", ylabel="Dist d", title="$gt Difference")
    hmd = heatmap!(ax_d, 0:nc-1, 0:nd-1, ds, colorrange=(-dmax, dmax), colormap=:RdBu)
    Colorbar(fig1b[gi, 5], hmd, label="Δ")
end
save_figure(fig1b, "smallworld_q1_heatmaps.png")


#==============================================================================#
#   QUESTIONS 2-4: Parameter Transfer vs Proxy                                 #
#==============================================================================#

println("\n" * "=" ^ 80)
println("QUESTIONS 2-4: Parameter Transfer vs Proxy")
println("=" ^ 80)
println("Source n=$N_SMALL, Target n=$N_LARGE, $NUM_INSTANCES instances\n")

source_instances = Dict(gt => generate_instances(gt, N_SMALL, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 1)) for gt in graph_types)
target_instances = Dict(gt => generate_instances(gt, N_LARGE, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 2)) for gt in graph_types)
target_optimal = Dict(gt => [maxcut_optimal(inst.costs) for inst in target_instances[gt]]
    for gt in graph_types)

println("--- Edge and Optimal Statistics ---")
for gt in graph_types
    se = [inst.num_edges for inst in source_instances[gt]]
    te = [inst.num_edges for inst in target_instances[gt]]
    @printf("  %-5s  src_edges: %.1f±%.1f  tgt_edges: %.1f±%.1f  c*: %.1f±%.1f\n",
        gt, mean(se), std(se), mean(te), std(te), mean(target_optimal[gt]), std(target_optimal[gt]))
end

# Compute optimal QAOA on a few target instances for reference
println("\n--- Computing Optimal QAOA on Target (first 5 instances, p=1) ---")
for gt in graph_types
    opt_ratios = map(1:min(5, NUM_INSTANCES)) do i
        inst = target_instances[gt][i]
        rng_opt = MersenneTwister(SEED + 9999 + i)
        _, _, best_exp = optimize_qaoa_cd(inst.costs, N_LARGE, 1, 30; rng=rng_opt)
        best_exp / target_optimal[gt][i]
    end
    @printf("  %-5s  optimal_p1_ratio: %.4f±%.4f\n", gt, mean(opt_ratios), std(opt_ratios))
end

# Random baseline
println("\n--- Random Parameter Baseline ---")
for gt in graph_types
    rng_rand = MersenneTwister(SEED + 7777)
    rand_ratios = map(enumerate(target_instances[gt])) do (i, inst)
        # Average over 100 random parameter sets
        best = -Inf
        for _ in 1:100
            γ = [rand(rng_rand) * 1.6]
            β = [rand(rng_rand) * π/2]
            exp_val = qaoa_expectation(inst.costs, N_LARGE, γ, β)
            best = max(best, exp_val)
        end
        best / target_optimal[gt][i]
    end
    print_stats("$gt random (best of 100, p=1)", rand_ratios)
end

# Methods: Transfer_p1 (coord descent), Transfer_ramp (for p>1 fair comparison), PaperProxy
# For p=1: Transfer uses coord descent on source, Proxy uses 2D grid
# For p>1: Both Transfer and Proxy use the same linear ramp grid, but Transfer
#           evaluates on real QAOA and Proxy evaluates on the proxy.

methods_p1 = ["Transfer", "PaperProxy"]
methods_ramp = ["Transfer (ramp)", "PaperProxy (ramp)"]

results = Dict(gt => Dict{String, Dict{Int, Vector{Float64}}}() for gt in graph_types)
for gt in graph_types
    for m in vcat(methods_p1, methods_ramp)
        results[gt][m] = Dict{Int, Vector{Float64}}()
    end
end

for gt in graph_types
    println("\n--- $gt ---")

    for p in P_VALUES
        println("  p=$p:")

        if p == 1
            # === Transfer (coord descent on source, median params on target) ===
            rng_opt = MersenneTwister(SEED + hash(gt) + p * 100)
            src_params = map(source_instances[gt]) do inst
                γs, βs, _ = optimize_qaoa_cd(inst.costs, N_SMALL, 1, N_RESTARTS; rng=rng_opt)
                (γs=γs, βs=βs)
            end
            med_γ = [median([sp.γs[1] for sp in src_params])]
            med_β = [median([sp.βs[1] for sp in src_params])]

            transfer_ratios = map(enumerate(target_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_LARGE, med_γ, med_β) / target_optimal[gt][i]
            end
            results[gt]["Transfer"][1] = transfer_ratios
            print_stats("Transfer (coord descent)", transfer_ratios)

            # === PaperProxy (2D grid) ===
            proxy_ratios = map(enumerate(target_instances[gt])) do (i, inst)
                proxy = make_proxy(gt, inst.num_edges, N_LARGE)
                homodist = cpu_compute_homodist(proxy)
                P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, 1)
                qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_optimal[gt][i]
            end
            results[gt]["PaperProxy"][1] = proxy_ratios
            print_stats("PaperProxy (2D grid)", proxy_ratios)

        else  # p > 1: fair comparison with same linear ramp grid
            # === Transfer (linear ramp, real QAOA grid search on source, median on target) ===
            src_ramp_params = map(source_instances[gt]) do inst
                ramp_params, _, _, _ = optimize_qaoa_linramp(inst.costs, N_SMALL, p)
                ramp_params  # (γ₁, γ_f, β₁, β_f) tuple
            end

            # Median of ramp parameters
            med_γ₁ = median([sp[1] for sp in src_ramp_params])
            med_γf = median([sp[2] for sp in src_ramp_params])
            med_β₁ = median([sp[3] for sp in src_ramp_params])
            med_βf = median([sp[4] for sp in src_ramp_params])
            γs_t, βs_t = linear_ramp(med_γ₁, med_γf, med_β₁, med_βf, p)

            ramp_transfer_ratios = map(enumerate(target_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_LARGE, γs_t .* π, βs_t .* π) / target_optimal[gt][i]
            end
            results[gt]["Transfer (ramp)"][p] = ramp_transfer_ratios
            print_stats("Transfer (linear ramp)", ramp_transfer_ratios)

            # === PaperProxy (same linear ramp grid, proxy evaluation) ===
            proxy_ramp_ratios = map(enumerate(target_instances[gt])) do (i, inst)
                proxy = make_proxy(gt, inst.num_edges, N_LARGE)
                homodist = cpu_compute_homodist(proxy)
                P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
                qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_optimal[gt][i]
            end
            results[gt]["PaperProxy (ramp)"][p] = proxy_ramp_ratios
            print_stats("PaperProxy (linear ramp)", proxy_ramp_ratios)
        end
    end
end


#==============================================================================#
#   CROSS-TYPE TRANSFER                                                        #
#==============================================================================#

println("\n" * "=" ^ 80)
println("CROSS-TYPE PARAMETER TRANSFER (p=1)")
println("=" ^ 80)

cross_results = Dict{String, Dict{String, Vector{Float64}}}()
for src in graph_types
    cross_results[src] = Dict{String, Vector{Float64}}()
    rng_c = MersenneTwister(SEED + hash(src) + 100)
    src_params = map(source_instances[src]) do inst
        γs, βs, _ = optimize_qaoa_cd(inst.costs, N_SMALL, 1, N_RESTARTS; rng=rng_c)
        (γs=γs, βs=βs)
    end
    med_γ = [median([sp.γs[1] for sp in src_params])]
    med_β = [median([sp.βs[1] for sp in src_params])]

    for tgt in graph_types
        ratios = map(enumerate(target_instances[tgt])) do (i, inst)
            qaoa_expectation(inst.costs, N_LARGE, med_γ, med_β) / target_optimal[tgt][i]
        end
        cross_results[src][tgt] = ratios
        print_stats("  $(src)→$(tgt)", ratios)
    end
end


#==============================================================================#
#   PLOTS                                                                      #
#==============================================================================#

println("\n--- Final Plots ---")

# Main comparison
fig2 = Figure(size=(400 * length(graph_types), 500))
method_colors = Dict(
    "Transfer" => :steelblue, "PaperProxy" => :coral,
    "Transfer (ramp)" => :steelblue, "PaperProxy (ramp)" => :coral)

for (gi, gt) in enumerate(graph_types)
    ax = Axis(fig2[1, gi],
        xlabel="QAOA Depth p", ylabel="Approx Ratio",
        title="$gt (n=$N_SMALL→$N_LARGE)", xticks=P_VALUES)

    for p in P_VALUES
        meths = p == 1 ? methods_p1 : methods_ramp
        nm = length(meths)
        tw = 0.6; sw = tw / nm

        for (mi, method) in enumerate(meths)
            offset = (mi - (nm + 1) / 2) * sw
            vals = results[gt][method][p]
            boxplot!(ax, fill(Float64(p) + offset, length(vals)), vals,
                color=method_colors[method], width=sw * 0.8,
                label=(p == P_VALUES[1] ? method : nothing))
        end
    end
    if gi == 1; axislegend(ax, position=:rb); end
end

Label(fig2[0, :],
    "Transfer vs PaperProxy(p_eff)\n$NUM_INSTANCES instances, same linear ramp grid for p>1",
    fontsize=14, font=:bold)
save_figure(fig2, "smallworld_q2_q3_q4_comparison.png")

# Cross-type heatmap (p=1 only)
fig3 = Figure(size=(500, 400))
ax3 = Axis(fig3[1, 1],
    xlabel="Target", ylabel="Source",
    title="Cross-Type Transfer (p=1, n=$N_SMALL→$N_LARGE)",
    xticks=(1:3, graph_types), yticks=(1:3, graph_types))

mat = [mean(cross_results[s][t]) for s in graph_types, t in graph_types]
hm = heatmap!(ax3, 1:3, 1:3, mat, colormap=:YlOrRd, colorrange=(minimum(mat) - 0.02, maximum(mat) + 0.02))
for i in 1:3, j in 1:3
    text!(ax3, j, i, text=@sprintf("%.3f", mat[i,j]),
        align=(:center, :center), fontsize=14, color=:black)
end
Colorbar(fig3[1, 2], hm, label="Approx Ratio")
save_figure(fig3, "smallworld_cross_type_transfer.png")


#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY TABLE")
println("=" ^ 80)

println("\np=1 (Transfer: coord descent on source; Proxy: 2D grid):")
@printf("%-8s  %-15s  %-15s\n", "Graph", "Transfer", "PaperProxy")
println("-" ^ 42)
for gt in graph_types
    t = results[gt]["Transfer"][1]
    p_r = results[gt]["PaperProxy"][1]
    @printf("%-8s  %.4f±%.4f   %.4f±%.4f\n", gt, mean(t), std(t), mean(p_r), std(p_r))
end

for p in P_VALUES
    p == 1 && continue
    println("\np=$p (Both use linear ramp grid, $GRID_SIZE_RAMP^4=$(GRID_SIZE_RAMP^4) points):")
    @printf("%-8s  %-18s  %-18s\n", "Graph", "Transfer(ramp)", "PaperProxy(ramp)")
    println("-" ^ 50)
    for gt in graph_types
        t = results[gt]["Transfer (ramp)"][p]
        p_r = results[gt]["PaperProxy (ramp)"][p]
        @printf("%-8s  %.4f±%.4f      %.4f±%.4f\n", gt, mean(t), std(t), mean(p_r), std(p_r))
    end
end

println("\nCross-type transfer (p=1):")
@printf("%-8s", "Src\\Tgt")
for gt in graph_types; @printf("  %-10s", gt); end
println()
for src in graph_types
    @printf("%-8s", src)
    for tgt in graph_types
        @printf("  %.4f    ", mean(cross_results[src][tgt]))
    end
    println()
end

println("\nDone!")
