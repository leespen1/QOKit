#=
fitted_triangle_investigation.jl — Fit TriangleProxy to BA/WS empirical
N(c';d,c) and evaluate QAOA performance vs Transfer and PaperProxy.

Research question: Can TriangleProxy, fitted to empirical homodist averaged
over many graph instances, outperform PaperProxy(p_eff) on non-ER graphs?

Approach:
  1. Generate 50 instances of each graph type at n=10
  2. Compute averaged empirical N(c';d,c) for each type
  3. Fit TriangleProxy parameters via grid search + refinement
  4. Evaluate: Transfer vs PaperProxy(p_eff) vs FittedTriangleProxy at p=1,2,3

Started: 2026-04-01
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

# Fitting configuration
const N_FIT = 10                  # n for fitting homodist
const NUM_FIT_INSTANCES = 50      # instances for averaging homodist

# Evaluation configuration
const N_SMALL = 9                 # source graph size
const N_LARGE = 12                # target graph size
const NUM_INSTANCES = 20          # instances per graph type for evaluation

# Graph parameters
const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

# QAOA parameters
const P_VALUES = [1, 2, 3]
const N_RESTARTS = 15
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10

const SEED = 42

#==============================================================================#
#                    FITTING HELPERS                                            #
#==============================================================================#

"""
Generate a grid of IntuitiveTriangleProxy instances for fitting.
Returns (proxies, param_matrix) where param_matrix is N×4.
"""
function generate_proxy_grid(num_constraints, num_qubits;
    height_range=range(0.2, 4.0, length=8),
    center_range=range(0.3, 0.7, length=6),
    left_angle_range=range(0.05, 0.45, length=6),
    right_angle_range=range(0.05, 0.45, length=6))

    proxies = TriangleProxy[]
    params = Vector{Float64}[]

    for h in height_range, c in center_range, la in left_angle_range, ra in right_angle_range
        try
            proxy = IntuitiveTriangleProxy(num_constraints, num_qubits, h, c, la, ra)
            push!(proxies, proxy)
            push!(params, [h, c, la, ra])
        catch
            continue  # skip invalid parameter combinations
        end
    end

    return proxies, reduce(hcat, params)'
end

"""
Local refinement of TriangleProxy parameters around a starting point.
Uses shrinking random perturbations (same strategy as sendai_opt).
"""
function refine_proxy_params(num_constraints, num_qubits, init_params, target_homodist;
    max_iter=2000, shrink_after=5, stop_after=30, init_sd=0.15)

    best_params = copy(init_params)
    best_proxy = IntuitiveTriangleProxy(num_constraints, num_qubits, best_params...)
    best_mse = cpu_multi_proxy_mse([best_proxy], target_homodist; normalize=true)[1]

    bounds = [
        (0.05, 8.0),   # height_adjustment
        (0.1, 0.9),    # center_adjustment
        (0.01, 0.49),  # left_angle
        (0.01, 0.49),  # right_angle
    ]

    sd = [init_sd * (b[2] - b[1]) for b in bounds]
    consecutive_fails = 0
    last_helpful_delta = zeros(4)

    for iter in 1:max_iter
        # Generate perturbation
        delta = randn(4) .* sd
        # Reuse helpful direction 30% of the time
        if !all(last_helpful_delta .== 0) && rand() < 0.3
            delta = last_helpful_delta .* (0.5 + rand())
        end

        trial_params = best_params .+ delta
        # Clamp to bounds
        for i in 1:4
            trial_params[i] = clamp(trial_params[i], bounds[i]...)
        end

        try
            trial_proxy = IntuitiveTriangleProxy(num_constraints, num_qubits, trial_params...)
            trial_mse = cpu_multi_proxy_mse([trial_proxy], target_homodist; normalize=true)[1]

            if trial_mse < best_mse
                last_helpful_delta = trial_params .- best_params
                best_params = trial_params
                best_mse = trial_mse
                consecutive_fails = 0
            else
                consecutive_fails += 1
            end
        catch
            consecutive_fails += 1
        end

        if consecutive_fails >= stop_after
            break
        elseif consecutive_fails >= shrink_after && consecutive_fails % shrink_after == 0
            sd .*= 0.7
        end
    end

    return best_params, best_mse
end


#==============================================================================#
#                    QAOA OPTIMIZATION HELPERS                                  #
#==============================================================================#

"""Optimize QAOA by coordinate descent (for p=1)."""
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

"""Optimize via proxy grid search (p=1: 2D grid, p>1: linear ramp grid)."""
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

function make_proxy(graph_type::String, num_edges::Int, num_vertices::Int)
    p_e = graph_type == "ER" ? ER_P_EDGE : effective_edge_probability(num_edges, num_vertices)
    PaperProxy(num_edges, num_vertices, p_e)
end

function print_stats(label, values)
    @printf("  %-45s  mean=%.4f  std=%.4f  [%.4f, %.4f]\n",
        label, mean(values), std(values), minimum(values), maximum(values))
end


#==============================================================================#
#   STEP 1: COMPUTE EMPIRICAL HOMODIST FOR EACH GRAPH TYPE                     #
#==============================================================================#

println("=" ^ 80)
println("STEP 1: Compute Empirical Homodist (n=$N_FIT, $NUM_FIT_INSTANCES instances)")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]

# Generate fitting instances
fit_instances = Dict(gt => generate_instances(gt, N_FIT, NUM_FIT_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 10)) for gt in graph_types)

# Find max edges for consistent homodist shape
max_edges_fit = maximum(inst.num_edges for gt in graph_types for inst in fit_instances[gt])
println("Max edges across all types: $max_edges_fit")

# Compute averaged empirical homodist
empirical_homodists = Dict{String, Array{Float64, 3}}()
for gt in graph_types
    println("  Computing homodist for $gt...")
    homodists = map(fit_instances[gt]) do inst
        get_homogeneous_distribution_from_costs_direct(
            inst.costs, inst.num_edges, inst.num_vertices;
            max_num_edges=max_edges_fit)
    end
    empirical_homodists[gt] = average_distributions(homodists)
    avg_edges = round(Int, mean(inst.num_edges for inst in fit_instances[gt]))
    @printf("    %s: avg_edges=%d, homodist shape=%s\n", gt, avg_edges, size(empirical_homodists[gt]))
end


#==============================================================================#
#   STEP 2: FIT TRIANGLEPROXY TO EACH GRAPH TYPE                               #
#==============================================================================#

println("\n" * "=" ^ 80)
println("STEP 2: Fit TriangleProxy Parameters")
println("=" ^ 80)

fitted_params = Dict{String, Vector{Float64}}()
fitted_mse = Dict{String, Float64}()

for gt in graph_types
    println("\n--- Fitting TriangleProxy for $gt ---")
    avg_edges = round(Int, mean(inst.num_edges for inst in fit_instances[gt]))
    target = empirical_homodists[gt]

    # Phase 1: Grid search
    println("  Phase 1: Grid search...")
    proxies, param_matrix = generate_proxy_grid(avg_edges, N_FIT)
    println("  Generated $(length(proxies)) proxy configurations")

    # Pad target to match proxy homodist size (proxy uses avg_edges, target uses max_edges_fit)
    proxy_ref = cpu_compute_homodist(proxies[1])
    target_padded, _ = pad_to_match(target, proxy_ref)
    # Now use proxy_ref shape for target
    target_for_mse = zeros(size(proxy_ref))
    sz = min.(size(target), size(target_for_mse))
    target_for_mse[1:sz[1], 1:sz[2], 1:sz[3]] = target[1:sz[1], 1:sz[2], 1:sz[3]]

    mse_values = cpu_multi_proxy_mse(proxies, target_for_mse; normalize=true)
    best_grid_idx = argmin(mse_values)
    best_grid_params = param_matrix[best_grid_idx, :]
    best_grid_mse = mse_values[best_grid_idx]
    @printf("  Best grid:  h=%.3f c=%.3f la=%.3f ra=%.3f  MSE=%.6e\n",
        best_grid_params..., best_grid_mse)

    # Phase 2: Local refinement
    println("  Phase 2: Local refinement...")
    refined_params, refined_mse = refine_proxy_params(
        avg_edges, N_FIT, collect(best_grid_params), target_for_mse;
        max_iter=3000, shrink_after=5, stop_after=40)
    @printf("  Refined:    h=%.3f c=%.3f la=%.3f ra=%.3f  MSE=%.6e\n",
        refined_params..., refined_mse)

    # Also compare with default TriangleProxy
    default_proxy = OldTriangleProxy(avg_edges, N_FIT)
    default_mse = cpu_multi_proxy_mse([default_proxy], target_for_mse; normalize=true)[1]
    @printf("  Default:    MSE=%.6e  (%.1fx worse than fitted)\n",
        default_mse, default_mse / refined_mse)

    # Compare with PaperProxy
    p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(avg_edges, N_FIT)
    paper_proxy = PaperProxy(avg_edges, N_FIT, p_eff)
    paper_homodist = cpu_compute_homodist(paper_proxy)
    # Pad paper homodist to match target shape
    paper_for_mse = zeros(size(target_for_mse))
    sz_p = min.(size(paper_homodist), size(paper_for_mse))
    paper_for_mse[1:sz_p[1], 1:sz_p[2], 1:sz_p[3]] = paper_homodist[1:sz_p[1], 1:sz_p[2], 1:sz_p[3]]
    # Normalize and compute MSE
    paper_for_mse_norm = paper_for_mse ./ max(sum(paper_for_mse), 1e-10)
    target_for_mse_norm = target_for_mse ./ max(sum(target_for_mse), 1e-10)
    paper_mse_val = sum((paper_for_mse_norm .- target_for_mse_norm).^2) / length(target_for_mse_norm)
    @printf("  PaperProxy: MSE=%.6e  (%.1fx vs fitted)\n",
        paper_mse_val, paper_mse_val / refined_mse)

    fitted_params[gt] = refined_params
    fitted_mse[gt] = refined_mse
end


#==============================================================================#
#   STEP 3: EVALUATE QAOA PERFORMANCE                                          #
#==============================================================================#

println("\n" * "=" ^ 80)
println("STEP 3: Approximation Ratio Comparison")
println("=" ^ 80)
println("Source n=$N_SMALL, Target n=$N_LARGE, $NUM_INSTANCES instances per type\n")

# Generate evaluation instances
source_instances = Dict(gt => generate_instances(gt, N_SMALL, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 1)) for gt in graph_types)
target_instances = Dict(gt => generate_instances(gt, N_LARGE, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 2)) for gt in graph_types)
target_optimal = Dict(gt => [maxcut_optimal(inst.costs) for inst in target_instances[gt]]
    for gt in graph_types)

println("--- Target Statistics ---")
for gt in graph_types
    te = [inst.num_edges for inst in target_instances[gt]]
    @printf("  %-5s  tgt_edges: %.1f±%.1f  c*: %.1f±%.1f\n",
        gt, mean(te), std(te), mean(target_optimal[gt]), std(target_optimal[gt]))
end

# Store results: results[gt][method][p] = Vector{Float64}
results = Dict(gt => Dict{String, Dict{Int, Vector{Float64}}}() for gt in graph_types)
methods = ["Transfer", "PaperProxy", "FittedTriangle"]
for gt in graph_types, m in methods
    results[gt][m] = Dict{Int, Vector{Float64}}()
end

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
            # Linear ramp grid search on source
            src_ramp_params = map(source_instances[gt]) do inst
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
            med_ramp = Tuple(median([sp[i] for sp in src_ramp_params]) for i in 1:4)
            γs_t, βs_t = linear_ramp(med_ramp..., p)

            transfer_ratios = map(enumerate(target_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_LARGE, γs_t .* π, βs_t .* π) / target_optimal[gt][i]
            end
        end
        results[gt]["Transfer"][p] = transfer_ratios
        print_stats("Transfer", transfer_ratios)

        # === PAPERPROXY (p_eff) ===
        proxy_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            proxy = make_proxy(gt, inst.num_edges, N_LARGE)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["PaperProxy"][p] = proxy_ratios
        print_stats("PaperProxy(p_eff)", proxy_ratios)

        # === FITTED TRIANGLEPROXY ===
        triangle_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            # Create fitted proxy for this target's edge count
            fp = fitted_params[gt]
            proxy = IntuitiveTriangleProxy(inst.num_edges, N_LARGE, fp...)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["FittedTriangle"][p] = triangle_ratios
        print_stats("FittedTriangle", triangle_ratios)
    end
end


#==============================================================================#
#   STEP 4: PEARSON CORRELATION COMPARISON                                     #
#==============================================================================#

println("\n" * "=" ^ 80)
println("STEP 4: Pearson Correlation (Fitted vs PaperProxy vs Empirical)")
println("=" ^ 80)

for gt in graph_types
    println("\n--- $gt ---")
    avg_edges = round(Int, mean(inst.num_edges for inst in fit_instances[gt]))
    target = empirical_homodists[gt]

    # PaperProxy homodist
    p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(avg_edges, N_FIT)
    paper_proxy = PaperProxy(avg_edges, N_FIT, p_eff)
    paper_hd = cpu_compute_homodist(paper_proxy)

    # Fitted TriangleProxy homodist
    fp = fitted_params[gt]
    fitted_proxy = IntuitiveTriangleProxy(avg_edges, N_FIT, fp...)
    fitted_hd = cpu_compute_homodist(fitted_proxy)

    # Compute P(c') for weighting
    num_bs = 1 << N_FIT
    P_vals = zeros(max_edges_fit + 1)
    for inst in fit_instances[gt]
        for x in 1:num_bs
            c = Int(inst.costs[x])
            if c <= max_edges_fit
                P_vals[c + 1] += 1.0 / (num_bs * NUM_FIT_INSTANCES)
            end
        end
    end

    for (label, hd) in [("PaperProxy", paper_hd), ("FittedTriangle", fitted_hd)]
        ph, eh = pad_to_match(hd, target)
        corrs = get_pearson_correlation_coefficients(ph, eh)
        num_costs = min(length(corrs), length(P_vals))
        valid = [i for i in 1:num_costs if !isnan(corrs[i]) && P_vals[i] > 0.01]
        if !isempty(valid)
            weighted_corr = sum(corrs[i] * P_vals[i] for i in valid) / sum(P_vals[i] for i in valid)
            min_corr = minimum(corrs[i] for i in valid)
            @printf("  %-20s  P-weighted_corr=%.4f  min_corr(dominant)=%.4f  n_dominant=%d\n",
                label, weighted_corr, min_corr, length(valid))
        end
    end
end


#==============================================================================#
#   PLOTS                                                                      #
#==============================================================================#

println("\n--- Plotting ---")

# Main comparison: box plots
method_colors = Dict(
    "Transfer" => :steelblue,
    "PaperProxy" => :coral,
    "FittedTriangle" => :mediumseagreen)

fig = Figure(size=(400 * length(graph_types), 550))

for (gi, gt) in enumerate(graph_types)
    ax = Axis(fig[1, gi],
        xlabel="QAOA Depth p", ylabel="Approx Ratio",
        title="$gt (n=$N_SMALL→$N_LARGE)", xticks=P_VALUES)

    for p in P_VALUES
        nm = length(methods)
        tw = 0.7; sw = tw / nm
        for (mi, method) in enumerate(methods)
            offset = (mi - (nm + 1) / 2) * sw
            vals = results[gt][method][p]
            boxplot!(ax, fill(Float64(p) + offset, length(vals)), vals,
                color=method_colors[method], width=sw * 0.8,
                label=(p == P_VALUES[1] ? method : nothing))
        end
    end
    if gi == 1; axislegend(ax, position=:rb); end
end

Label(fig[0, :],
    "Transfer vs PaperProxy(p_eff) vs FittedTriangleProxy\n$NUM_INSTANCES instances, fitted on $NUM_FIT_INSTANCES instances at n=$N_FIT",
    fontsize=14, font=:bold)
save_figure(fig, "fitted_triangle_comparison.png")

# Homodist comparison heatmaps (for the middle cost slice)
fig2 = Figure(size=(1400, 350 * length(graph_types)))
for (gi, gt) in enumerate(graph_types)
    avg_edges = round(Int, mean(inst.num_edges for inst in fit_instances[gt]))
    target = empirical_homodists[gt]
    c_prime = size(target, 1) ÷ 2
    emp_s = target[c_prime, :, :]
    nd, nc = size(emp_s)
    vmax = max(maximum(emp_s), 1e-10)

    # Empirical
    ax_e = Axis(fig2[gi, 1], xlabel="Cost c", ylabel="Dist d",
        title="$gt Empirical, c'=$(c_prime-1)")
    heatmap!(ax_e, 0:nc-1, 0:nd-1, emp_s, colorrange=(0, vmax), colormap=:viridis)

    # PaperProxy
    p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(avg_edges, N_FIT)
    paper_hd = cpu_compute_homodist(PaperProxy(avg_edges, N_FIT, p_eff))
    ph_padded, _ = pad_to_match(paper_hd, target)
    ps = ph_padded[c_prime, 1:nd, 1:nc]
    ax_p = Axis(fig2[gi, 2], xlabel="Cost c", ylabel="Dist d",
        title="$gt PaperProxy(p_eff)")
    heatmap!(ax_p, 0:nc-1, 0:nd-1, ps, colorrange=(0, vmax), colormap=:viridis)

    # FittedTriangle
    fp = fitted_params[gt]
    fitted_hd = cpu_compute_homodist(IntuitiveTriangleProxy(avg_edges, N_FIT, fp...))
    fh_padded, _ = pad_to_match(fitted_hd, target)
    fs = fh_padded[c_prime, 1:nd, 1:nc]
    ax_f = Axis(fig2[gi, 3], xlabel="Cost c", ylabel="Dist d",
        title="$gt FittedTriangle")
    heatmap!(ax_f, 0:nc-1, 0:nd-1, fs, colorrange=(0, vmax), colormap=:viridis)

    # Difference: Fitted - Empirical
    ds = fs .- emp_s
    dmax = max(maximum(abs.(ds)), 1e-10)
    ax_d = Axis(fig2[gi, 4], xlabel="Cost c", ylabel="Dist d",
        title="$gt Fitted−Empirical")
    hmd = heatmap!(ax_d, 0:nc-1, 0:nd-1, ds, colorrange=(-dmax, dmax), colormap=:RdBu)
    Colorbar(fig2[gi, 5], hmd, label="Δ")
end
save_figure(fig2, "fitted_triangle_heatmaps.png")


#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

println("\nFitted TriangleProxy Parameters:")
@printf("%-8s  %-10s  %-10s  %-10s  %-10s  %-12s\n",
    "Graph", "height", "center", "left_a", "right_a", "MSE")
println("-" ^ 65)
for gt in graph_types
    fp = fitted_params[gt]
    @printf("%-8s  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %.6e\n",
        gt, fp..., fitted_mse[gt])
end

println("\nApproximation Ratios (mean ± std):")
for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-8s  %-20s  %-20s  %-20s\n", "Graph", "Transfer", "PaperProxy", "FittedTriangle")
    println("  " * "-" ^ 72)
    for gt in graph_types
        t = results[gt]["Transfer"][p]
        pp = results[gt]["PaperProxy"][p]
        ft = results[gt]["FittedTriangle"][p]
        @printf("  %-8s  %.4f±%.4f        %.4f±%.4f        %.4f±%.4f\n",
            gt, mean(t), std(t), mean(pp), std(pp), mean(ft), std(ft))
    end
end

# Win/loss summary
println("\nWin Counts (FittedTriangle vs PaperProxy, per instance):")
for gt in ["BA", "WS"]
    for p in P_VALUES
        ft_wins = sum(results[gt]["FittedTriangle"][p] .> results[gt]["PaperProxy"][p])
        pp_wins = sum(results[gt]["PaperProxy"][p] .> results[gt]["FittedTriangle"][p])
        ties = NUM_INSTANCES - ft_wins - pp_wins
        @printf("  %-5s p=%d: FittedTriangle=%d  PaperProxy=%d  ties=%d\n",
            gt, p, ft_wins, pp_wins, ties)
    end
end

println("\nDone!")
