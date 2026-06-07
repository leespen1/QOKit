#=
gaussian_proxy_investigation.jl — Test Gaussian blob proxy on non-ER graphs.

The Gaussian proxy represents N(c';d,c) as a Gaussian shape in (d,c) space
for each c' slice. More expressive than TriangleProxy (which failed at p>1)
but still parametric.

Approach:
1. Compute empirical homodist for fitting instances
2. Fit GaussianProxy parameters by grid search to minimize MSE
3. Use the fitted proxy with consistent Gaussian P(c') for QAOA optimization
4. Compare against Transfer, PaperProxy, SampN+EmpP

Key question: does a richer functional form survive multi-layer propagation?

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_TARGET = 14
const N_SOURCE = 9
const NUM_FIT_INSTANCES = 10       # instances for homodist fitting
const NUM_EVAL_INSTANCES = 10      # instances for QAOA evaluation
const SAMPLES_PER_COST = 20
const P_VALUES = [1, 3, 5]
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const SEED = 42

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

#==============================================================================#
#                    PROXY FITTING                                              #
#==============================================================================#

"""
Fit GaussianProxy parameters to empirical homodist via grid search + random perturbation.
"""
function fit_gaussian_proxy(empirical_homodist::Array{Float64, 3}, num_edges::Int, num_qubits::Int)
    m = num_edges
    n = num_qubits

    function objective(params)
        center_target, sigma_base, sigma_scale, height_base, height_scale, height_power, center_bias = params

        proxy = GaussianProxy(m, n;
            center_target=clamp(center_target, 0.1, 0.9),
            sigma_base=clamp(sigma_base, 0.1, 10.0),
            sigma_scale=clamp(sigma_scale, 0.0, 5.0),
            height_base=clamp(height_base, 0.01, 10.0),
            height_scale=clamp(height_scale, 0.0, 5.0),
            height_power=clamp(height_power, 0.1, 2.0),
            center_bias=clamp(center_bias, -1.0, 1.0))

        proxy_homodist = cpu_compute_homodist(proxy)

        # Trim to match dimensions
        s = min.(size(proxy_homodist), size(empirical_homodist))
        p_view = @view proxy_homodist[1:s[1], 1:s[2], 1:s[3]]
        e_view = @view empirical_homodist[1:s[1], 1:s[2], 1:s[3]]

        # Normalize for fair comparison
        p_sum = max(sum(p_view), 1e-10)
        e_sum = max(sum(e_view), 1e-10)

        mse = 0.0
        for i in eachindex(p_view)
            diff = p_view[i] / p_sum - e_view[i] / e_sum
            mse += diff * diff
        end
        return mse / length(p_view)
    end

    # Bounds for each parameter
    lo = [0.1, 0.1, 0.0, 0.01, 0.0, 0.1, -1.0]
    hi = [0.9, 10.0, 5.0, 10.0, 5.0, 2.0, 1.0]

    # Phase 1: Coarse grid search
    best_params = [0.5, 1.0, 1.0, 1.0, 1.0, 0.5, 0.0]
    best_mse = objective(best_params)

    rng = MersenneTwister(SEED + 999)
    for _ in 1:200
        params = [lo[i] + rand(rng) * (hi[i] - lo[i]) for i in 1:7]
        mse = objective(params)
        if mse < best_mse
            best_mse = mse
            best_params = params
        end
    end

    # Phase 2: Local refinement via random perturbation (smart random search)
    step_scale = 0.2
    consecutive_failures = 0
    direction = zeros(7)

    for iter in 1:500
        # Random perturbation, reuse successful direction
        perturbation = randn(rng, 7) .* step_scale .* (hi .- lo)
        if consecutive_failures < 3 && any(direction .!= 0)
            perturbation .+= 0.3 .* direction
        end

        candidate = best_params .+ perturbation
        candidate = clamp.(candidate, lo, hi)
        mse = objective(candidate)

        if mse < best_mse
            direction = candidate .- best_params
            best_mse = mse
            best_params = candidate
            consecutive_failures = 0
        else
            consecutive_failures += 1
            if consecutive_failures > 10
                step_scale *= 0.8
                consecutive_failures = 0
            end
        end
    end

    proxy = GaussianProxy(m, n;
        center_target=clamp(best_params[1], 0.1, 0.9),
        sigma_base=clamp(best_params[2], 0.1, 10.0),
        sigma_scale=clamp(best_params[3], 0.0, 5.0),
        height_base=clamp(best_params[4], 0.01, 10.0),
        height_scale=clamp(best_params[5], 0.0, 5.0),
        height_power=clamp(best_params[6], 0.1, 2.0),
        center_bias=clamp(best_params[7], -1.0, 1.0))

    return proxy, best_mse
end

#==============================================================================#
#                    SHARED HELPERS                                             #
#==============================================================================#

function generate_instances(graph_type::String, n::Int, num::Int; rng=Random.default_rng())
    map(1:num) do _
        if graph_type == "ER"; generate_er_instance(n, ER_P_EDGE; rng)
        elseif graph_type == "BA"; generate_ba_instance(n, BA_M_ATTACH; rng)
        elseif graph_type == "WS"; generate_ws_instance(n, WS_K, WS_P_REWIRE; rng)
        else; error("Unknown: $graph_type"); end
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
println("Gaussian Blob Proxy Investigation")
println("n=$N_TARGET, p=$P_VALUES")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]
methods = ["Random", "Transfer", "PaperProxy", "SampN+EmpP", "GaussianFit+GaussP", "GaussianFit+EmpP"]

all_results = Dict{String, Dict{String, Dict{Int, Vector{Float64}}}}()
for gt in graph_types
    all_results[gt] = Dict(m => Dict{Int, Vector{Float64}}() for m in methods)
end

for gt in graph_types
    println("\n" * "=" ^ 40)
    println("  $gt (n=$N_TARGET)")
    println("=" ^ 40)

    # Generate instances
    eval_rng = MersenneTwister(SEED + hash(gt) + N_TARGET * 100)
    eval_insts = generate_instances(gt, N_TARGET, NUM_EVAL_INSTANCES; rng=eval_rng)
    eval_optimal = [maxcut_optimal(inst.costs) for inst in eval_insts]

    fit_rng = MersenneTwister(SEED + hash(gt) + N_TARGET * 200)
    fit_insts = generate_instances(gt, N_TARGET, NUM_FIT_INSTANCES; rng=fit_rng)

    src_rng = MersenneTwister(SEED + hash(gt) + 888)
    source_insts = generate_instances(gt, N_SOURCE, NUM_EVAL_INSTANCES; rng=src_rng)

    # Compute empirical homodist for fitting
    max_edges = maximum(inst.num_edges for inst in fit_insts)
    println("  Computing sampled homodist for fitting (max_edges=$max_edges)...")
    t_hd = @elapsed begin
        sampled_homodists = map(fit_insts) do inst
            get_homogeneous_distribution_from_costs_sampled(
                inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
                max_num_edges=max_edges,
                rng=MersenneTwister(SEED + hash(gt) + hash(inst.num_edges) + N_TARGET))
        end
        empirical_avg = average_distributions(sampled_homodists)
    end
    println("  Homodist computed in $(round(t_hd, digits=2))s")

    # Fit Gaussian proxy
    println("  Fitting GaussianProxy...")
    m_fit = max_edges
    t_fit = @elapsed begin
        fitted_proxy, fit_mse = fit_gaussian_proxy(empirical_avg, m_fit, N_TARGET)
    end
    @printf("  Fitted in %.2fs, MSE=%.6f\n", t_fit, fit_mse)
    @printf("  Params: center=%.3f σ_base=%.3f σ_scale=%.3f h_base=%.3f h_scale=%.3f h_pow=%.3f bias=%.3f\n",
        fitted_proxy.center_target, fitted_proxy.sigma_base, fitted_proxy.sigma_scale,
        fitted_proxy.height_base, fitted_proxy.height_scale, fitted_proxy.height_power,
        fitted_proxy.center_bias)

    # Compute proxy homodist and P values
    gauss_homodist = cpu_compute_homodist(fitted_proxy)
    gauss_P = [P_cost_distribution(fitted_proxy, c) for c in 0:fitted_proxy.num_constraints]

    # Empirical P
    P_emp = compute_empirical_P(fit_insts, max_edges)
    P_emp_trimmed = P_emp[1:min(m_fit+1, length(P_emp))]
    if length(P_emp_trimmed) < m_fit + 1
        P_emp_trimmed = vcat(P_emp_trimmed, zeros(m_fit + 1 - length(P_emp_trimmed)))
    end

    for p in P_VALUES
        println("\n  p=$p:")

        # === RANDOM ===
        random_ratios = [mean(inst.costs) / eval_optimal[i] for (i, inst) in enumerate(eval_insts)]
        all_results[gt]["Random"][p] = random_ratios
        print_stats("Random baseline", random_ratios)

        # === TRANSFER ===
        t_γ, t_β = transfer_optimize(source_insts, N_SOURCE, p)
        transfer_ratios = [qaoa_expectation(inst.costs, N_TARGET, t_γ, t_β) / eval_optimal[i]
                           for (i, inst) in enumerate(eval_insts)]
        all_results[gt]["Transfer"][p] = transfer_ratios
        print_stats("Transfer", transfer_ratios)

        # === PAPERPROXY ===
        pp_ratios = map(enumerate(eval_insts)) do (i, inst)
            p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, N_TARGET)
            proxy = PaperProxy(inst.num_edges, N_TARGET, p_eff)
            hd = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy(hd, P_vals, N_TARGET, p)
            qaoa_expectation(inst.costs, N_TARGET, best_γ, best_β) / eval_optimal[i]
        end
        all_results[gt]["PaperProxy"][p] = pp_ratios
        print_stats("PaperProxy", pp_ratios)

        # === SAMPN+EMPP ===
        best_γ_s, best_β_s, _ = optimize_via_proxy(empirical_avg, P_emp_trimmed, N_TARGET, p)
        samp_ratios = [qaoa_expectation(inst.costs, N_TARGET, best_γ_s, best_β_s) / eval_optimal[i]
                       for (i, inst) in enumerate(eval_insts)]
        all_results[gt]["SampN+EmpP"][p] = samp_ratios
        print_stats("SampN+EmpP", samp_ratios)

        # === GAUSSIAN + GAUSSIAN P (consistent) ===
        # Trim homodist to match proxy dimensions
        s_gauss = min.(size(gauss_homodist), size(empirical_avg))
        gauss_hd_trimmed = gauss_homodist[1:s_gauss[1], 1:s_gauss[2], 1:s_gauss[3]]
        gauss_P_trimmed = gauss_P[1:min(length(gauss_P), s_gauss[1])]
        if length(gauss_P_trimmed) < s_gauss[1]
            gauss_P_trimmed = vcat(gauss_P_trimmed, zeros(s_gauss[1] - length(gauss_P_trimmed)))
        end

        best_γ_g, best_β_g, _ = optimize_via_proxy(gauss_hd_trimmed, gauss_P_trimmed, N_TARGET, p)
        gauss_ratios = [qaoa_expectation(inst.costs, N_TARGET, best_γ_g, best_β_g) / eval_optimal[i]
                       for (i, inst) in enumerate(eval_insts)]
        all_results[gt]["GaussianFit+GaussP"][p] = gauss_ratios
        print_stats("GaussianFit+GaussP", gauss_ratios)

        # === GAUSSIAN + EMPIRICAL P (test consistency) ===
        best_γ_ge, best_β_ge, _ = optimize_via_proxy(gauss_hd_trimmed, P_emp_trimmed[1:s_gauss[1]], N_TARGET, p)
        gauss_emp_ratios = [qaoa_expectation(inst.costs, N_TARGET, best_γ_ge, best_β_ge) / eval_optimal[i]
                       for (i, inst) in enumerate(eval_insts)]
        all_results[gt]["GaussianFit+EmpP"][p] = gauss_emp_ratios
        print_stats("GaussianFit+EmpP", gauss_emp_ratios)
    end
end

#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY: Approximation Ratios at n=$N_TARGET")
println("=" ^ 80)

for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-5s", "Graph")
    for m in methods; @printf("  %-18s", m); end
    println()
    println("  " * "-" ^ (7 + 20 * length(methods)))
    for gt in graph_types
        @printf("  %-5s", gt)
        for m in methods
            if haskey(all_results[gt][m], p)
                vals = all_results[gt][m][p]
                @printf("  %.4f±%.4f    ", mean(vals), std(vals))
            else
                @printf("  %-18s", "—")
            end
        end
        println()
    end
end

println("\nDone!")
