#=
hybrid_warmstart_investigation.jl — Test hybrid proxy warmstart + local
refinement on real QAOA.

The idea: proxy gives a cheap-to-evaluate approximation of the QAOA landscape.
Use it to find a good starting point, then refine that starting point with a
few coordinate descent steps on the REAL QAOA. This combines:
- Proxy's cheap global search
- Real QAOA's accuracy for local refinement

Methods compared:
  1. Transfer (baseline — optimize source, median to target)
  2. PaperProxy only (proxy grid search, no refinement)
  3. PaperProxy + Refine (proxy warmstart + coord descent on target)
  4. Random + Refine (random start + same coord descent — control)

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

const P_VALUES = [1, 2, 3, 5]
const N_RESTARTS = 15
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10

# Refinement parameters (how many coord descent steps on real QAOA)
const REFINE_RESTARTS = 3        # restarts around proxy optimum
const REFINE_PERTURBATION = 0.2  # perturbation radius for restarts
const RANDOM_RESTARTS = 10       # restarts for random baseline

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

"""Coord descent starting from a specific initial point (for warmstart)."""
function refine_qaoa_from(costs, n, γs_init, βs_init)
    current_γs = copy(γs_init)
    current_βs = copy(βs_init)
    p = length(γs_init)
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
    return current_γs, current_βs, current_exp
end

"""Warmstart: proxy optimum + perturbation restarts, refined on real QAOA."""
function hybrid_optimize(costs, n, proxy_γs, proxy_βs; n_restarts=REFINE_RESTARTS,
    perturb=REFINE_PERTURBATION, rng=Random.default_rng())
    p = length(proxy_γs)
    best_γs, best_βs, best_exp = refine_qaoa_from(costs, n, proxy_γs, proxy_βs)

    for _ in 1:n_restarts
        # Perturbed restart around proxy optimum
        trial_γs = proxy_γs .+ randn(rng, p) .* perturb
        trial_βs = proxy_βs .+ randn(rng, p) .* perturb
        trial_γs = max.(trial_γs, 0)
        trial_βs = clamp.(trial_βs, 0, π/2)
        ref_γs, ref_βs, ref_exp = refine_qaoa_from(costs, n, trial_γs, trial_βs)
        if ref_exp > best_exp
            best_exp = ref_exp
            best_γs, best_βs = ref_γs, ref_βs
        end
    end
    return best_γs, best_βs, best_exp
end

function optimize_via_proxy(homodist, P_vals, n, p)
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
        K = gs^4
        γ_matrix = zeros(K, p)
        β_matrix = zeros(K, p)
        idx = 0
        for γ₁ in range(0.02, 0.40, length=gs),
            γ_f in range(0.10, 0.70, length=gs),
            β₁ in range(0.05, 0.45, length=gs),
            β_f in range(0.01, 0.25, length=gs)
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
println("HYBRID WARMSTART: Proxy + Local Refinement")
println("Source n=$N_SMALL, Target n=$N_LARGE, $NUM_INSTANCES instances")
println("Depths: $(P_VALUES)")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]
methods = ["Transfer", "PaperProxy", "Proxy+Refine", "Random+Refine"]

source_instances = Dict(gt => generate_instances(gt, N_SMALL, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 1)) for gt in graph_types)
target_instances = Dict(gt => generate_instances(gt, N_LARGE, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 2)) for gt in graph_types)
target_optimal = Dict(gt => [maxcut_optimal(inst.costs) for inst in target_instances[gt]]
    for gt in graph_types)

results = Dict(gt => Dict(m => Dict{Int, Vector{Float64}}() for m in methods) for gt in graph_types)

for gt in graph_types
    println("\n--- $gt ---")

    for p in P_VALUES
        println("  p=$p:")

        # === TRANSFER ===
        if p == 1
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
        else
            src_ramp = map(source_instances[gt]) do inst
                best_params = (0.0, 0.0, 0.0, 0.0); best_exp = -Inf
                gs = GRID_SIZE_RAMP
                for γ₁ in range(0.02, 0.40, length=gs),
                    γ_f in range(0.10, 0.70, length=gs),
                    β₁ in range(0.05, 0.45, length=gs),
                    β_f in range(0.01, 0.25, length=gs)
                    γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, p)
                    exp_val = qaoa_expectation(inst.costs, N_SMALL, γs_pi .* π, βs_pi .* π)
                    if exp_val > best_exp; best_exp = exp_val; best_params = (γ₁, γ_f, β₁, β_f); end
                end
                best_params
            end
            med_ramp = Tuple(median([sp[i] for sp in src_ramp]) for i in 1:4)
            γs_t, βs_t = linear_ramp(med_ramp..., p)
            transfer_ratios = map(enumerate(target_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_LARGE, γs_t .* π, βs_t .* π) / target_optimal[gt][i]
            end
        end
        results[gt]["Transfer"][p] = transfer_ratios
        @printf("    Transfer:      mean=%.4f  std=%.4f\n", mean(transfer_ratios), std(transfer_ratios))

        # === PAPERPROXY ONLY ===
        proxy_only_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, N_LARGE)
            proxy = PaperProxy(inst.num_edges, N_LARGE, p_eff)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["PaperProxy"][p] = proxy_only_ratios
        @printf("    PaperProxy:    mean=%.4f  std=%.4f\n", mean(proxy_only_ratios), std(proxy_only_ratios))

        # === PROXY + REFINE (hybrid) ===
        hybrid_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, N_LARGE)
            proxy = PaperProxy(inst.num_edges, N_LARGE, p_eff)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            proxy_γ, proxy_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
            # Refine on real QAOA
            rng_h = MersenneTwister(SEED + i * 1000)
            ref_γ, ref_β, ref_exp = hybrid_optimize(inst.costs, N_LARGE, proxy_γ, proxy_β; rng=rng_h)
            ref_exp / target_optimal[gt][i]
        end
        results[gt]["Proxy+Refine"][p] = hybrid_ratios
        @printf("    Proxy+Refine:  mean=%.4f  std=%.4f\n", mean(hybrid_ratios), std(hybrid_ratios))

        # === RANDOM + REFINE (control) ===
        random_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            rng_r = MersenneTwister(SEED + i * 2000)
            best_exp = -Inf
            best_γs = zeros(p)
            best_βs = zeros(p)
            for _ in 1:RANDOM_RESTARTS
                init_γ = rand(rng_r, p) .* 1.6
                init_β = rand(rng_r, p) .* (π/2)
                ref_γ, ref_β, ref_exp = refine_qaoa_from(inst.costs, N_LARGE, init_γ, init_β)
                if ref_exp > best_exp
                    best_exp = ref_exp
                    best_γs, best_βs = ref_γ, ref_β
                end
            end
            best_exp / target_optimal[gt][i]
        end
        results[gt]["Random+Refine"][p] = random_ratios
        @printf("    Random+Refine: mean=%.4f  std=%.4f\n", mean(random_ratios), std(random_ratios))
    end
end


#==============================================================================#
#   PLOTS                                                                      #
#==============================================================================#

println("\n--- Plotting ---")

method_colors = Dict(
    "Transfer" => :steelblue,
    "PaperProxy" => :coral,
    "Proxy+Refine" => :mediumseagreen,
    "Random+Refine" => :gray60)

fig = Figure(size=(500 * length(graph_types), 550))
for (gi, gt) in enumerate(graph_types)
    ax = Axis(fig[1, gi],
        xlabel="QAOA Depth p", ylabel="Approx Ratio",
        title="$gt (n=$N_SMALL→$N_LARGE)", xticks=P_VALUES)

    for p in P_VALUES
        nm = length(methods); tw = 0.8; sw = tw / nm
        for (mi, method) in enumerate(methods)
            offset = (mi - (nm + 1) / 2) * sw
            vals = results[gt][method][p]
            boxplot!(ax, fill(Float64(p) + offset, length(vals)), vals,
                color=method_colors[method], width=sw * 0.8,
                label=(p == P_VALUES[1] ? method : nothing))
        end
    end
    if gi == 1; axislegend(ax, position=:rb, labelsize=9); end
end

Label(fig[0, :],
    "Hybrid Warmstart: Proxy + Local Refinement on Real QAOA\n$NUM_INSTANCES instances, $(REFINE_RESTARTS)+1 refinement runs, $RANDOM_RESTARTS random restarts",
    fontsize=14, font=:bold)
save_figure(fig, "hybrid_warmstart.png")


#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-8s", "Graph")
    for m in methods; @printf("  %-16s", m); end
    println()
    println("  " * "-" ^ (8 + 18 * length(methods)))
    for gt in graph_types
        @printf("  %-8s", gt)
        for m in methods
            vals = results[gt][m][p]
            @printf("  %.4f±%.4f  ", mean(vals), std(vals))
        end
        println()
    end
end

println("\n\nKey Comparison (Proxy+Refine vs Transfer, positive = Hybrid better):")
for p in P_VALUES
    @printf("  p=%d:", p)
    for gt in graph_types
        diff = mean(results[gt]["Proxy+Refine"][p]) - mean(results[gt]["Transfer"][p])
        @printf("  %s=%+.4f", gt, diff)
    end
    println()
end

println("\nProxy+Refine Improvement over PaperProxy:")
for p in P_VALUES
    @printf("  p=%d:", p)
    for gt in graph_types
        diff = mean(results[gt]["Proxy+Refine"][p]) - mean(results[gt]["PaperProxy"][p])
        @printf("  %s=%+.4f", gt, diff)
    end
    println()
end

println("\nDone!")
