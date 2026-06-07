#=
depth_scaling_investigation.jl — Depth scaling (p=1 to p=8) with linear ramp
for all graph types. How does approximation ratio grow with depth?

Tests Transfer vs PaperProxy(p_eff) at depths 1-8 using linear ramp schedules.
Characterizes how the proxy-vs-transfer gap changes as a function of depth.

Started: 2026-04-01
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

const P_VALUES = collect(1:8)
const N_RESTARTS = 15
const GRID_SIZE_RAMP = 10  # 10^4 = 10000 points per depth
const SEED = 42

#==============================================================================#
#                    HELPERS                                                    #
#==============================================================================#

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

function optimize_linramp_real(costs, n, p; grid_size=GRID_SIZE_RAMP)
    best_params = (0.0, 0.0, 0.0, 0.0)
    best_exp = -Inf
    for γ₁ in range(0.02, 0.40, length=grid_size),
        γ_f in range(0.10, 0.70, length=grid_size),
        β₁ in range(0.05, 0.45, length=grid_size),
        β_f in range(0.01, 0.25, length=grid_size)
        γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, p)
        exp_val = qaoa_expectation(costs, n, γs_pi .* π, βs_pi .* π)
        if exp_val > best_exp; best_exp = exp_val; best_params = (γ₁, γ_f, β₁, β_f); end
    end
    return best_params, best_exp
end

function optimize_via_proxy_ramp(homodist, P_vals, n, p; grid_size=GRID_SIZE_RAMP)
    K = grid_size^4
    γ_matrix = zeros(K, p)
    β_matrix = zeros(K, p)
    idx = 0
    for γ₁ in range(0.02, 0.40, length=grid_size),
        γ_f in range(0.10, 0.70, length=grid_size),
        β₁ in range(0.05, 0.45, length=grid_size),
        β_f in range(0.01, 0.25, length=grid_size)
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


#==============================================================================#
#   MAIN EXPERIMENT                                                            #
#==============================================================================#

println("=" ^ 80)
println("DEPTH SCALING: Transfer vs PaperProxy(p_eff)")
println("Source n=$N_SMALL, Target n=$N_LARGE, $NUM_INSTANCES instances")
println("Depths: $(P_VALUES)")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]

source_instances = Dict(gt => generate_instances(gt, N_SMALL, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 1)) for gt in graph_types)
target_instances = Dict(gt => generate_instances(gt, N_LARGE, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 2)) for gt in graph_types)
target_optimal = Dict(gt => [maxcut_optimal(inst.costs) for inst in target_instances[gt]]
    for gt in graph_types)

methods = ["Transfer", "PaperProxy"]
results = Dict(gt => Dict(m => Dict{Int, Vector{Float64}}() for m in methods) for gt in graph_types)

for gt in graph_types
    println("\n--- $gt ---")

    for p in P_VALUES
        println("  p=$p:")

        # === TRANSFER ===
        if p == 1
            rng_opt = MersenneTwister(SEED + hash(gt) + 100)
            src_params = map(source_instances[gt]) do inst
                γs, βs, _ = optimize_qaoa_cd(inst.costs, N_SMALL, 1, N_RESTARTS; rng=rng_opt)
                (γs=γs, βs=βs)
            end
            med_γ = [median([sp.γs[1] for sp in src_params])]
            med_β = [median([sp.βs[1] for sp in src_params])]
            transfer_ratios = map(enumerate(target_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_LARGE, med_γ, med_β) / target_optimal[gt][i]
            end
        else
            src_ramp = map(source_instances[gt]) do inst
                ramp_params, _ = optimize_linramp_real(inst.costs, N_SMALL, p)
                ramp_params
            end
            med_ramp = Tuple(median([sp[i] for sp in src_ramp]) for i in 1:4)
            γs_t, βs_t = linear_ramp(med_ramp..., p)
            transfer_ratios = map(enumerate(target_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_LARGE, γs_t .* π, βs_t .* π) / target_optimal[gt][i]
            end
        end
        results[gt]["Transfer"][p] = transfer_ratios
        @printf("    Transfer:   mean=%.4f  std=%.4f\n", mean(transfer_ratios), std(transfer_ratios))

        # === PAPERPROXY ===
        proxy_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, N_LARGE)
            proxy = PaperProxy(inst.num_edges, N_LARGE, p_eff)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy_ramp(homodist, P_vals, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["PaperProxy"][p] = proxy_ratios
        @printf("    PaperProxy: mean=%.4f  std=%.4f\n", mean(proxy_ratios), std(proxy_ratios))
    end
end


#==============================================================================#
#   PLOTS                                                                      #
#==============================================================================#

println("\n--- Plotting ---")

method_colors = Dict("Transfer" => :steelblue, "PaperProxy" => :coral)

fig = Figure(size=(500 * length(graph_types), 500))

for (gi, gt) in enumerate(graph_types)
    ax = Axis(fig[1, gi],
        xlabel="QAOA Depth p", ylabel="Approx Ratio",
        title="$gt (n=$N_SMALL→$N_LARGE)", xticks=P_VALUES)

    for method in methods
        means = [mean(results[gt][method][p]) for p in P_VALUES]
        stds = [std(results[gt][method][p]) for p in P_VALUES]
        band!(ax, P_VALUES, means .- stds, means .+ stds,
            color=(method_colors[method], 0.2))
        lines!(ax, P_VALUES, means, color=method_colors[method],
            linewidth=2, label=method)
        scatter!(ax, P_VALUES, means, color=method_colors[method], markersize=8)
    end
    if gi == 1; axislegend(ax, position=:rb); end
end

Label(fig[0, :],
    "Depth Scaling: Transfer vs PaperProxy(p_eff)\n$NUM_INSTANCES instances, linear ramp grid ($GRID_SIZE_RAMP^4 points)",
    fontsize=14, font=:bold)
save_figure(fig, "depth_scaling.png")

# Gap plot
fig2 = Figure(size=(600, 400))
ax2 = Axis(fig2[1, 1],
    xlabel="QAOA Depth p", ylabel="Transfer − PaperProxy Gap",
    title="Proxy-Transfer Gap vs Depth")

gt_colors = Dict("ER" => :steelblue, "BA" => :coral, "WS" => :mediumseagreen)
for gt in graph_types
    gaps = [mean(results[gt]["Transfer"][p]) - mean(results[gt]["PaperProxy"][p]) for p in P_VALUES]
    lines!(ax2, P_VALUES, gaps, color=gt_colors[gt], linewidth=2, label=gt)
    scatter!(ax2, P_VALUES, gaps, color=gt_colors[gt], markersize=8)
end
hlines!(ax2, [0], color=:gray, linestyle=:dash)
axislegend(ax2, position=:lt)
save_figure(fig2, "depth_scaling_gap.png")


#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

@printf("\n%-8s", "p")
for gt in graph_types
    @printf("  %-12s %-12s %-8s", "$gt Xfer", "$gt Proxy", "Gap")
end
println()
println("-" ^ (8 + 35 * length(graph_types)))

for p in P_VALUES
    @printf("%-8d", p)
    for gt in graph_types
        t_m = mean(results[gt]["Transfer"][p])
        p_m = mean(results[gt]["PaperProxy"][p])
        @printf("  %-12.4f %-12.4f %-8.4f", t_m, p_m, t_m - p_m)
    end
    println()
end

println("\nDone!")
