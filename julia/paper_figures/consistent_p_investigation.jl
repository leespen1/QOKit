#=
consistent_p_investigation.jl — Test whether pairing FittedTriangleProxy's
N(c';d,c) with a consistent P(c') resolves the depth-compounding failure.

Hypothesis: The fitted triangle's poor QAOA performance at depth>1 is caused
by inconsistency between its N(c';d,c) and its crude triangular P(c').
If we substitute a better P(c'), the fitted proxy might work.

Methods tested:
  1. Transfer (baseline)
  2. PaperProxy (N + P both from paper formula)
  3. FittedTriangle + TriangleP (original — N from fit, P from triangle formula)
  4. FittedTriangle + BinomialP (N from fit, P from PaperProxy's binomial)
  5. FittedTriangle + EmpiricalP (N from fit, P from empirical distribution)
  6. EmpiricalN + EmpiricalP (N from averaged instances, P from same instances)

Started: 2026-04-01
=#

include("smallworld_common.jl")
using Printf
using Distributions: Binomial, pdf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_FIT = 10
const NUM_FIT_INSTANCES = 50
const N_SMALL = 9
const N_LARGE = 12
const NUM_INSTANCES = 20

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

const P_VALUES = [1, 2, 3]
const N_RESTARTS = 15
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
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

function print_stats(label, values)
    @printf("  %-55s  mean=%.4f  std=%.4f\n", label, mean(values), std(values))
end

"""Compute empirical P(c') from a list of instances."""
function compute_empirical_P(instances, max_edges)
    P = zeros(max_edges + 1)
    n = instances[1].num_vertices
    num_bs = 1 << n
    for inst in instances
        for x in 1:num_bs
            c = Int(inst.costs[x])
            if c <= max_edges
                P[c + 1] += 1.0 / (num_bs * length(instances))
            end
        end
    end
    return P
end

"""Compute binomial P(c') for given edge count and probability."""
function compute_binomial_P(num_edges, p_edge)
    d = Binomial(num_edges, p_edge)
    return [pdf(d, c) for c in 0:num_edges]
end


#==============================================================================#
#   STEP 1: COMPUTE FITTED TRIANGLE PARAMS AND EMPIRICAL DATA                  #
#==============================================================================#

println("=" ^ 80)
println("STEP 1: Prepare Proxies and Distributions")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]

# Generate fitting instances
fit_instances = Dict(gt => generate_instances(gt, N_FIT, NUM_FIT_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 10)) for gt in graph_types)

max_edges_fit = maximum(inst.num_edges for gt in graph_types for inst in fit_instances[gt])

# Compute empirical homodist and P(c') for fitting
empirical_homodists_fit = Dict{String, Array{Float64, 3}}()
empirical_P_fit = Dict{String, Vector{Float64}}()

for gt in graph_types
    println("  Computing homodist for $gt...")
    homodists = map(fit_instances[gt]) do inst
        get_homogeneous_distribution_from_costs_direct(
            inst.costs, inst.num_edges, inst.num_vertices;
            max_num_edges=max_edges_fit)
    end
    empirical_homodists_fit[gt] = average_distributions(homodists)
    empirical_P_fit[gt] = compute_empirical_P(fit_instances[gt], max_edges_fit)
end

# Fit TriangleProxy (reuse parameters from fitted_triangle_investigation)
# These were found by grid search + refinement on n=10, 50 instances
fitted_params = Dict{String, Vector{Float64}}()

# Quick re-fit using same approach
for gt in graph_types
    println("  Fitting TriangleProxy for $gt...")
    avg_edges = round(Int, mean(inst.num_edges for inst in fit_instances[gt]))
    target = empirical_homodists_fit[gt]

    # Generate grid
    proxies = TriangleProxy[]
    params_list = Vector{Float64}[]
    for h in range(0.2, 6.0, length=10),
        c in range(0.3, 0.7, length=7),
        la in range(0.05, 0.45, length=7),
        ra in range(0.05, 0.45, length=7)
        try
            proxy = IntuitiveTriangleProxy(avg_edges, N_FIT, h, c, la, ra)
            push!(proxies, proxy)
            push!(params_list, [h, c, la, ra])
        catch; continue; end
    end

    # Pad target for MSE computation
    proxy_ref = cpu_compute_homodist(proxies[1])
    target_for_mse = zeros(size(proxy_ref))
    sz = min.(size(target), size(target_for_mse))
    target_for_mse[1:sz[1], 1:sz[2], 1:sz[3]] = target[1:sz[1], 1:sz[2], 1:sz[3]]

    mse_values = cpu_multi_proxy_mse(proxies, target_for_mse; normalize=true)
    best_idx = argmin(mse_values)
    best_params = params_list[best_idx]

    # Refine
    bounds = [(0.05, 8.0), (0.1, 0.9), (0.01, 0.49), (0.01, 0.49)]
    sd = [0.15 * (b[2] - b[1]) for b in bounds]
    consecutive_fails = 0
    best_mse = mse_values[best_idx]
    for _ in 1:2000
        delta = randn(4) .* sd
        trial = clamp.(best_params .+ delta, [b[1] for b in bounds], [b[2] for b in bounds])
        try
            proxy = IntuitiveTriangleProxy(avg_edges, N_FIT, trial...)
            mse = cpu_multi_proxy_mse([proxy], target_for_mse; normalize=true)[1]
            if mse < best_mse
                best_params = trial; best_mse = mse; consecutive_fails = 0
            else
                consecutive_fails += 1
            end
        catch; consecutive_fails += 1; end
        consecutive_fails >= 40 && break
        consecutive_fails >= 5 && consecutive_fails % 5 == 0 && (sd .*= 0.7)
    end

    fitted_params[gt] = best_params
    @printf("  %s: h=%.3f c=%.3f la=%.3f ra=%.3f  MSE=%.6e\n", gt, best_params..., best_mse)
end


#==============================================================================#
#   STEP 2: EVALUATE ALL METHODS                                               #
#==============================================================================#

println("\n" * "=" ^ 80)
println("STEP 2: Approximation Ratio Comparison")
println("=" ^ 80)

source_instances = Dict(gt => generate_instances(gt, N_SMALL, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 1)) for gt in graph_types)
target_instances = Dict(gt => generate_instances(gt, N_LARGE, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 2)) for gt in graph_types)
target_optimal = Dict(gt => [maxcut_optimal(inst.costs) for inst in target_instances[gt]]
    for gt in graph_types)

# Also compute empirical homodist and P for target-sized graphs (for EmpiricalN+EmpiricalP)
# Use separate instances to avoid data leakage
target_fit_instances = Dict(gt => generate_instances(gt, N_LARGE, NUM_FIT_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 50)) for gt in graph_types)

max_edges_target = maximum(inst.num_edges for gt in graph_types for inst in target_fit_instances[gt])

empirical_homodists_target = Dict{String, Array{Float64, 3}}()
empirical_P_target = Dict{String, Vector{Float64}}()

println("Computing empirical homodist at n=$N_LARGE...")
for gt in graph_types
    println("  $gt...")
    homodists = map(target_fit_instances[gt]) do inst
        get_homogeneous_distribution_from_costs_direct(
            inst.costs, inst.num_edges, inst.num_vertices;
            max_num_edges=max_edges_target)
    end
    empirical_homodists_target[gt] = average_distributions(homodists)
    empirical_P_target[gt] = compute_empirical_P(target_fit_instances[gt], max_edges_target)
end

methods = [
    "Transfer",
    "PaperProxy",
    "Fitted+TriP",
    "Fitted+BinP",
    "Fitted+EmpP",
    "EmpN+EmpP",
]

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
        print_stats("Transfer", transfer_ratios)

        # === PAPERPROXY ===
        proxy_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, N_LARGE)
            proxy = PaperProxy(inst.num_edges, N_LARGE, p_eff)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["PaperProxy"][p] = proxy_ratios
        print_stats("PaperProxy", proxy_ratios)

        # === FITTED TRIANGLE + TRIANGLE P (original) ===
        fitted_tri_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            fp = fitted_params[gt]
            proxy = IntuitiveTriangleProxy(inst.num_edges, N_LARGE, fp...)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["Fitted+TriP"][p] = fitted_tri_ratios
        print_stats("Fitted+TriP", fitted_tri_ratios)

        # === FITTED TRIANGLE + BINOMIAL P ===
        fitted_bin_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            fp = fitted_params[gt]
            proxy = IntuitiveTriangleProxy(inst.num_edges, N_LARGE, fp...)
            homodist = cpu_compute_homodist(proxy)
            p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, N_LARGE)
            P_vals = compute_binomial_P(inst.num_edges, p_eff)
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["Fitted+BinP"][p] = fitted_bin_ratios
        print_stats("Fitted+BinP", fitted_bin_ratios)

        # === FITTED TRIANGLE + EMPIRICAL P ===
        fitted_emp_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            fp = fitted_params[gt]
            proxy = IntuitiveTriangleProxy(inst.num_edges, N_LARGE, fp...)
            homodist = cpu_compute_homodist(proxy)
            # Use empirical P from target-sized fitting instances, truncated to inst's edge count
            P_full = empirical_P_target[gt]
            P_vals = P_full[1:min(inst.num_edges+1, length(P_full))]
            # Pad if needed
            if length(P_vals) < inst.num_edges + 1
                P_vals = vcat(P_vals, zeros(inst.num_edges + 1 - length(P_vals)))
            end
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["Fitted+EmpP"][p] = fitted_emp_ratios
        print_stats("Fitted+EmpP", fitted_emp_ratios)

        # === EMPIRICAL N + EMPIRICAL P ===
        emp_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            homodist = empirical_homodists_target[gt]
            P_full = empirical_P_target[gt]
            # Truncate P to match homodist size
            m_hd = size(homodist, 1) - 1
            P_vals = P_full[1:min(m_hd+1, length(P_full))]
            if length(P_vals) < m_hd + 1
                P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
            end
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["EmpN+EmpP"][p] = emp_ratios
        print_stats("EmpN+EmpP", emp_ratios)
    end
end


#==============================================================================#
#   PLOTS                                                                      #
#==============================================================================#

println("\n--- Plotting ---")

method_colors = Dict(
    "Transfer" => :steelblue,
    "PaperProxy" => :coral,
    "Fitted+TriP" => :mediumseagreen,
    "Fitted+BinP" => :orange,
    "Fitted+EmpP" => :purple,
    "EmpN+EmpP" => :gold,
)

fig = Figure(size=(500 * length(graph_types), 600))
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
    "N/P Consistency Test: Which P(c') pairs best with fitted N(c';d,c)?\n$NUM_INSTANCES eval instances, fitted on $NUM_FIT_INSTANCES instances at n=$N_FIT",
    fontsize=14, font=:bold)
save_figure(fig, "consistent_p_comparison.png")


#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-8s", "Graph")
    for m in methods; @printf("  %-14s", m); end
    println()
    println("  " * "-" ^ (8 + 16 * length(methods)))
    for gt in graph_types
        @printf("  %-8s", gt)
        for m in methods
            @printf("  %.4f±%.4f", mean(results[gt][m][p]), std(results[gt][m][p]))
        end
        println()
    end
end

# Key comparison: Fitted+BinP vs PaperProxy (is binomial P the fix?)
println("\n\nKey Question: Does Binomial P fix the fitted proxy?")
println("(Fitted+BinP - PaperProxy, positive = fitted wins)")
for p in P_VALUES
    println("  p=$p:")
    for gt in ["BA", "WS"]
        diff = mean(results[gt]["Fitted+BinP"][p]) - mean(results[gt]["PaperProxy"][p])
        @printf("    %s: %+.4f\n", gt, diff)
    end
end

println("\nDone!")
