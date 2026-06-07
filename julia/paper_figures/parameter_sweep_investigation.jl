#=
parameter_sweep_investigation.jl — Sweep BA/WS graph parameters to test
robustness of proxy and transfer performance.

Research question: Do the results from the smallworld investigation (BA m=2,
WS k=4 p=0.3) hold across the graph parameter space?

Sweeps:
  - BA: m_attach ∈ {1, 2, 3, 4}
  - WS: (k, p_rewire) ∈ {(2, 0.1), (4, 0.3), (6, 0.5)}

Compares Transfer vs PaperProxy(p_eff) at p=1,2,3.

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

const ER_P_EDGE = 0.5
const P_VALUES = [1, 2, 3]
const N_RESTARTS = 15
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const SEED = 42

# Parameter sweep configurations
const BA_CONFIGS = [1, 2, 3, 4]  # m_attach values
const WS_CONFIGS = [
    (k=2, p_r=0.1),
    (k=4, p_r=0.3),
    (k=6, p_r=0.5),
]

#==============================================================================#
#                    HELPERS (reused from smallworld)                           #
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

function print_stats(label, values)
    @printf("  %-50s  mean=%.4f  std=%.4f\n", label, mean(values), std(values))
end


#==============================================================================#
#   BA PARAMETER SWEEP                                                         #
#==============================================================================#

println("=" ^ 80)
println("BA PARAMETER SWEEP: m_attach ∈ $(BA_CONFIGS)")
println("=" ^ 80)
println("Source n=$N_SMALL, Target n=$N_LARGE, $NUM_INSTANCES instances\n")

ba_results = Dict{Int, Dict{String, Dict{Int, Vector{Float64}}}}()

for m_attach in BA_CONFIGS
    println("\n--- BA m_attach=$m_attach ---")
    ba_results[m_attach] = Dict("Transfer" => Dict{Int, Vector{Float64}}(),
                                "PaperProxy" => Dict{Int, Vector{Float64}}())

    rng_src = MersenneTwister(SEED + hash(m_attach) + 1)
    rng_tgt = MersenneTwister(SEED + hash(m_attach) + 2)

    source = [generate_ba_instance(N_SMALL, m_attach; rng=rng_src) for _ in 1:NUM_INSTANCES]
    target = [generate_ba_instance(N_LARGE, m_attach; rng=rng_tgt) for _ in 1:NUM_INSTANCES]
    target_opt = [maxcut_optimal(inst.costs) for inst in target]

    avg_src_edges = mean(inst.num_edges for inst in source)
    avg_tgt_edges = mean(inst.num_edges for inst in target)
    @printf("  src_edges: %.1f  tgt_edges: %.1f  c*: %.1f±%.1f\n",
        avg_src_edges, avg_tgt_edges, mean(target_opt), std(target_opt))

    for p in P_VALUES
        println("  p=$p:")

        # Transfer
        if p == 1
            rng_opt = MersenneTwister(SEED + hash(m_attach) + p * 100)
            src_params = map(source) do inst
                γs, βs, _ = optimize_qaoa_cd(inst.costs, N_SMALL, 1, N_RESTARTS; rng=rng_opt)
                (γs=γs, βs=βs)
            end
            med_γ = [median([sp.γs[1] for sp in src_params])]
            med_β = [median([sp.βs[1] for sp in src_params])]
            transfer_ratios = map(enumerate(target)) do (i, inst)
                qaoa_expectation(inst.costs, N_LARGE, med_γ, med_β) / target_opt[i]
            end
        else
            src_ramp = map(source) do inst
                best_params = (0.0, 0.0, 0.0, 0.0)
                best_exp = -Inf
                gs = GRID_SIZE_RAMP
                for γ₁ in range(0.02, 0.40, length=gs),
                    γ_f in range(0.10, 0.70, length=gs),
                    β₁ in range(0.05, 0.45, length=gs),
                    β_f in range(0.01, 0.25, length=gs)
                    γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, p)
                    exp_val = qaoa_expectation(inst.costs, N_SMALL, γs_pi .* π, βs_pi .* π)
                    if exp_val > best_exp
                        best_exp = exp_val
                        best_params = (γ₁, γ_f, β₁, β_f)
                    end
                end
                best_params
            end
            med_ramp = Tuple(median([sp[i] for sp in src_ramp]) for i in 1:4)
            γs_t, βs_t = linear_ramp(med_ramp..., p)
            transfer_ratios = map(enumerate(target)) do (i, inst)
                qaoa_expectation(inst.costs, N_LARGE, γs_t .* π, βs_t .* π) / target_opt[i]
            end
        end
        ba_results[m_attach]["Transfer"][p] = transfer_ratios
        print_stats("Transfer", transfer_ratios)

        # PaperProxy(p_eff)
        proxy_ratios = map(enumerate(target)) do (i, inst)
            p_eff = effective_edge_probability(inst.num_edges, N_LARGE)
            proxy = PaperProxy(inst.num_edges, N_LARGE, p_eff)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_opt[i]
        end
        ba_results[m_attach]["PaperProxy"][p] = proxy_ratios
        print_stats("PaperProxy(p_eff)", proxy_ratios)
    end
end


#==============================================================================#
#   WS PARAMETER SWEEP                                                         #
#==============================================================================#

println("\n" * "=" ^ 80)
println("WS PARAMETER SWEEP: (k, p_rewire) ∈ $(WS_CONFIGS)")
println("=" ^ 80)

ws_results = Dict{String, Dict{String, Dict{Int, Vector{Float64}}}}()

for cfg in WS_CONFIGS
    label = "k=$(cfg.k),p=$(cfg.p_r)"
    println("\n--- WS $label ---")
    ws_results[label] = Dict("Transfer" => Dict{Int, Vector{Float64}}(),
                             "PaperProxy" => Dict{Int, Vector{Float64}}())

    rng_src = MersenneTwister(SEED + hash(label) + 1)
    rng_tgt = MersenneTwister(SEED + hash(label) + 2)

    source = [generate_ws_instance(N_SMALL, cfg.k, cfg.p_r; rng=rng_src) for _ in 1:NUM_INSTANCES]
    target = [generate_ws_instance(N_LARGE, cfg.k, cfg.p_r; rng=rng_tgt) for _ in 1:NUM_INSTANCES]
    target_opt = [maxcut_optimal(inst.costs) for inst in target]

    avg_src_edges = mean(inst.num_edges for inst in source)
    avg_tgt_edges = mean(inst.num_edges for inst in target)
    @printf("  src_edges: %.1f  tgt_edges: %.1f  c*: %.1f±%.1f\n",
        avg_src_edges, avg_tgt_edges, mean(target_opt), std(target_opt))

    for p in P_VALUES
        println("  p=$p:")

        # Transfer
        if p == 1
            rng_opt = MersenneTwister(SEED + hash(label) + p * 100)
            src_params = map(source) do inst
                γs, βs, _ = optimize_qaoa_cd(inst.costs, N_SMALL, 1, N_RESTARTS; rng=rng_opt)
                (γs=γs, βs=βs)
            end
            med_γ = [median([sp.γs[1] for sp in src_params])]
            med_β = [median([sp.βs[1] for sp in src_params])]
            transfer_ratios = map(enumerate(target)) do (i, inst)
                qaoa_expectation(inst.costs, N_LARGE, med_γ, med_β) / target_opt[i]
            end
        else
            src_ramp = map(source) do inst
                best_params = (0.0, 0.0, 0.0, 0.0)
                best_exp = -Inf
                gs = GRID_SIZE_RAMP
                for γ₁ in range(0.02, 0.40, length=gs),
                    γ_f in range(0.10, 0.70, length=gs),
                    β₁ in range(0.05, 0.45, length=gs),
                    β_f in range(0.01, 0.25, length=gs)
                    γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, p)
                    exp_val = qaoa_expectation(inst.costs, N_SMALL, γs_pi .* π, βs_pi .* π)
                    if exp_val > best_exp
                        best_exp = exp_val
                        best_params = (γ₁, γ_f, β₁, β_f)
                    end
                end
                best_params
            end
            med_ramp = Tuple(median([sp[i] for sp in src_ramp]) for i in 1:4)
            γs_t, βs_t = linear_ramp(med_ramp..., p)
            transfer_ratios = map(enumerate(target)) do (i, inst)
                qaoa_expectation(inst.costs, N_LARGE, γs_t .* π, βs_t .* π) / target_opt[i]
            end
        end
        ws_results[label]["Transfer"][p] = transfer_ratios
        print_stats("Transfer", transfer_ratios)

        # PaperProxy(p_eff)
        proxy_ratios = map(enumerate(target)) do (i, inst)
            p_eff = effective_edge_probability(inst.num_edges, N_LARGE)
            proxy = PaperProxy(inst.num_edges, N_LARGE, p_eff)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_opt[i]
        end
        ws_results[label]["PaperProxy"][p] = proxy_ratios
        print_stats("PaperProxy(p_eff)", proxy_ratios)
    end
end


#==============================================================================#
#   PLOTS                                                                      #
#==============================================================================#

println("\n--- Plotting ---")

method_colors = Dict("Transfer" => :steelblue, "PaperProxy" => :coral)
methods_list = ["Transfer", "PaperProxy"]

# BA sweep plot
fig_ba = Figure(size=(400 * length(BA_CONFIGS), 500))
for (ci, m_attach) in enumerate(BA_CONFIGS)
    ax = Axis(fig_ba[1, ci],
        xlabel="QAOA Depth p", ylabel="Approx Ratio",
        title="BA m=$m_attach", xticks=P_VALUES)
    for p in P_VALUES
        nm = 2; tw = 0.6; sw = tw / nm
        for (mi, method) in enumerate(methods_list)
            offset = (mi - 1.5) * sw
            vals = ba_results[m_attach][method][p]
            boxplot!(ax, fill(Float64(p) + offset, length(vals)), vals,
                color=method_colors[method], width=sw * 0.8,
                label=(p == 1 ? method : nothing))
        end
    end
    if ci == 1; axislegend(ax, position=:rb); end
end
Label(fig_ba[0, :],
    "BA Parameter Sweep: Transfer vs PaperProxy(p_eff)\nn=$N_SMALL→$N_LARGE, $NUM_INSTANCES instances",
    fontsize=14, font=:bold)
save_figure(fig_ba, "parameter_sweep_ba.png")

# WS sweep plot
ws_labels = ["k=$(c.k),p=$(c.p_r)" for c in WS_CONFIGS]
fig_ws = Figure(size=(400 * length(WS_CONFIGS), 500))
for (ci, label) in enumerate(ws_labels)
    ax = Axis(fig_ws[1, ci],
        xlabel="QAOA Depth p", ylabel="Approx Ratio",
        title="WS $label", xticks=P_VALUES)
    for p in P_VALUES
        nm = 2; tw = 0.6; sw = tw / nm
        for (mi, method) in enumerate(methods_list)
            offset = (mi - 1.5) * sw
            vals = ws_results[label][method][p]
            boxplot!(ax, fill(Float64(p) + offset, length(vals)), vals,
                color=method_colors[method], width=sw * 0.8,
                label=(p == 1 ? method : nothing))
        end
    end
    if ci == 1; axislegend(ax, position=:rb); end
end
Label(fig_ws[0, :],
    "WS Parameter Sweep: Transfer vs PaperProxy(p_eff)\nn=$N_SMALL→$N_LARGE, $NUM_INSTANCES instances",
    fontsize=14, font=:bold)
save_figure(fig_ws, "parameter_sweep_ws.png")


#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY: BA PARAMETER SWEEP")
println("=" ^ 80)
for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-12s  %-20s  %-20s  %-8s\n", "m_attach", "Transfer", "PaperProxy", "Gap")
    println("  " * "-" ^ 65)
    for m_attach in BA_CONFIGS
        t = ba_results[m_attach]["Transfer"][p]
        pp = ba_results[m_attach]["PaperProxy"][p]
        gap = mean(t) - mean(pp)
        @printf("  %-12d  %.4f±%.4f        %.4f±%.4f        %+.4f\n",
            m_attach, mean(t), std(t), mean(pp), std(pp), gap)
    end
end

println("\n" * "=" ^ 80)
println("SUMMARY: WS PARAMETER SWEEP")
println("=" ^ 80)
for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-15s  %-20s  %-20s  %-8s\n", "(k, p_rewire)", "Transfer", "PaperProxy", "Gap")
    println("  " * "-" ^ 68)
    for label in ws_labels
        t = ws_results[label]["Transfer"][p]
        pp = ws_results[label]["PaperProxy"][p]
        gap = mean(t) - mean(pp)
        @printf("  %-15s  %.4f±%.4f        %.4f±%.4f        %+.4f\n",
            label, mean(t), std(t), mean(pp), std(pp), gap)
    end
end

# Density comparison
println("\n" * "=" ^ 80)
println("EDGE DENSITY (p_eff) FOR EACH CONFIG")
println("=" ^ 80)
max_edges_12 = N_LARGE * (N_LARGE - 1) ÷ 2
println("  max possible edges for n=$N_LARGE: $max_edges_12")
for m_attach in BA_CONFIGS
    rng_tmp = MersenneTwister(999)
    avg_edges = mean(generate_ba_instance(N_LARGE, m_attach; rng=rng_tmp).num_edges for _ in 1:10)
    @printf("  BA m=%d:  avg_edges=%.1f  p_eff=%.3f\n", m_attach, avg_edges, avg_edges / max_edges_12)
end
for cfg in WS_CONFIGS
    rng_tmp = MersenneTwister(999)
    avg_edges = mean(generate_ws_instance(N_LARGE, cfg.k, cfg.p_r; rng=rng_tmp).num_edges for _ in 1:10)
    @printf("  WS k=%d,p=%.1f:  avg_edges=%.1f  p_eff=%.3f\n",
        cfg.k, cfg.p_r, avg_edges, avg_edges / max_edges_12)
end

println("\nDone!")
