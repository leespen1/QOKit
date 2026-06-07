#=
systematic_proxy_evaluation.jl — Head-to-head comparison of TriangleProxy vs
NormalProxy on non-ER graphs with consistent N/P pairing.

Key design decisions:
  - Each proxy uses its OWN P(c') (consistency principle)
  - Proxies are fitted to multi-instance averaged homodist at fitting size
  - Evaluation on separate target instances
  - Compare: Transfer, PaperProxy, FittedTriangle+TriP, FittedNormal+NormP,
    EmpN+EmpP, and (new) SampledN+EmpP

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf
using Distributions: Binomial, pdf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_FIT = 10                    # size for fitting proxies
const NUM_FIT_INSTANCES = 50        # instances for fitting
const N_TARGET = 12                 # size for evaluation
const NUM_TARGET_INSTANCES = 20     # instances for evaluation
const N_SOURCE = 9                  # size for transfer source

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

const P_VALUES = [1, 2, 3, 5]
const N_RESTARTS = 15
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const SEED = 42

const SAMPLES_PER_COST = 10        # for sampled homodist

#==============================================================================#
#                    HELPERS                                                    #
#==============================================================================#

function optimize_qaoa_cd(costs, n, p, n_restarts; rng=Random.default_rng())
    best_exp = -Inf; best_γs = zeros(p); best_βs = zeros(p)
    for _ in 1:n_restarts
        γs = rand(rng, p) .* 1.6; βs = rand(rng, p) .* (π/2)
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

function print_stats(label, values)
    @printf("  %-55s  mean=%.4f  std=%.4f\n", label, mean(values), std(values))
end


#==============================================================================#
#   STEP 1: FIT PROXIES                                                        #
#==============================================================================#

println("=" ^ 80)
println("STEP 1: Fit TriangleProxy and NormalProxy to Multi-Instance Homodist")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]

fit_instances = Dict(gt => generate_instances(gt, N_FIT, NUM_FIT_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 10)) for gt in graph_types)

max_edges_fit = maximum(inst.num_edges for gt in graph_types for inst in fit_instances[gt])

# Compute empirical homodist for fitting
empirical_homodists_fit = Dict{String, Array{Float64, 3}}()
for gt in graph_types
    println("  Computing homodist for $gt at n=$N_FIT...")
    homodists = map(fit_instances[gt]) do inst
        get_homogeneous_distribution_from_costs_direct(
            inst.costs, inst.num_edges, inst.num_vertices;
            max_num_edges=max_edges_fit)
    end
    empirical_homodists_fit[gt] = average_distributions(homodists)
end

# Fit TriangleProxy
fitted_tri_params = Dict{String, Vector{Float64}}()
for gt in graph_types
    println("  Fitting TriangleProxy for $gt...")
    avg_edges = round(Int, mean(inst.num_edges for inst in fit_instances[gt]))
    target = empirical_homodists_fit[gt]

    proxies = TriangleProxy[]; params_list = Vector{Float64}[]
    for h in range(0.2, 6.0, length=10),
        c in range(0.3, 0.7, length=7),
        la in range(0.05, 0.45, length=7),
        ra in range(0.05, 0.45, length=7)
        try
            proxy = IntuitiveTriangleProxy(avg_edges, N_FIT, h, c, la, ra)
            push!(proxies, proxy); push!(params_list, [h, c, la, ra])
        catch; continue; end
    end

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
    consecutive_fails = 0; best_mse = mse_values[best_idx]
    for _ in 1:2000
        delta = randn(4) .* sd
        trial = clamp.(best_params .+ delta, [b[1] for b in bounds], [b[2] for b in bounds])
        try
            proxy = IntuitiveTriangleProxy(avg_edges, N_FIT, trial...)
            mse = cpu_multi_proxy_mse([proxy], target_for_mse; normalize=true)[1]
            if mse < best_mse; best_params = trial; best_mse = mse; consecutive_fails = 0
            else; consecutive_fails += 1; end
        catch; consecutive_fails += 1; end
        consecutive_fails >= 40 && break
        consecutive_fails >= 5 && consecutive_fails % 5 == 0 && (sd .*= 0.7)
    end
    fitted_tri_params[gt] = best_params
    @printf("  %s TriangleProxy: h=%.3f c=%.3f la=%.3f ra=%.3f  MSE=%.6e\n", gt, best_params..., best_mse)
end

# Fit NormalProxy
fitted_norm_params = Dict{String, Vector{Float64}}()
for gt in graph_types
    println("  Fitting NormalProxy for $gt...")
    avg_edges = round(Int, mean(inst.num_edges for inst in fit_instances[gt]))
    target = empirical_homodists_fit[gt]

    # Grid search over (cost_mean, cov_1, cov_2)
    best_params = [Float64(avg_edges)/2, 1.0, 1.0]
    best_mse = Inf

    for cm in range(avg_edges * 0.3, avg_edges * 0.7, length=10),
        c1 in range(0.1, 5.0, length=10),
        c2 in range(0.1, 5.0, length=10)
        try
            proxy = NormalProxy(avg_edges, N_FIT, cm, c1, c2)
            homodist = cpu_compute_homodist(proxy)
            sz = min.(size(target), size(homodist))
            target_crop = target[1:sz[1], 1:sz[2], 1:sz[3]]
            homodist_crop = homodist[1:sz[1], 1:sz[2], 1:sz[3]]
            # Normalize both for fair comparison
            tn = target_crop ./ max(sum(target_crop), 1e-10)
            hn = homodist_crop ./ max(sum(homodist_crop), 1e-10)
            mse = sum((tn .- hn).^2) / length(tn)
            if mse < best_mse; best_mse = mse; best_params = [cm, c1, c2]; end
        catch; continue; end
    end

    # Refine
    bounds = [(avg_edges * 0.1, avg_edges * 0.9), (0.01, 10.0), (0.01, 10.0)]
    sd = [0.1 * (b[2] - b[1]) for b in bounds]
    consecutive_fails = 0
    for _ in 1:2000
        delta = randn(3) .* sd
        trial = clamp.(best_params .+ delta, [b[1] for b in bounds], [b[2] for b in bounds])
        try
            proxy = NormalProxy(avg_edges, N_FIT, trial...)
            homodist = cpu_compute_homodist(proxy)
            sz = min.(size(target), size(homodist))
            target_crop = target[1:sz[1], 1:sz[2], 1:sz[3]]
            homodist_crop = homodist[1:sz[1], 1:sz[2], 1:sz[3]]
            tn = target_crop ./ max(sum(target_crop), 1e-10)
            hn = homodist_crop ./ max(sum(homodist_crop), 1e-10)
            mse = sum((tn .- hn).^2) / length(tn)
            if mse < best_mse; best_params = trial; best_mse = mse; consecutive_fails = 0
            else; consecutive_fails += 1; end
        catch; consecutive_fails += 1; end
        consecutive_fails >= 40 && break
        consecutive_fails >= 5 && consecutive_fails % 5 == 0 && (sd .*= 0.7)
    end
    fitted_norm_params[gt] = best_params
    @printf("  %s NormalProxy: cm=%.3f c1=%.3f c2=%.3f  MSE=%.6e\n", gt, best_params..., best_mse)
end


#==============================================================================#
#   STEP 2: EVALUATE ALL METHODS                                               #
#==============================================================================#

println("\n" * "=" ^ 80)
println("STEP 2: Approximation Ratio Comparison")
println("=" ^ 80)

source_instances = Dict(gt => generate_instances(gt, N_SOURCE, NUM_TARGET_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 1)) for gt in graph_types)
target_instances = Dict(gt => generate_instances(gt, N_TARGET, NUM_TARGET_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 2)) for gt in graph_types)
target_optimal = Dict(gt => [maxcut_optimal(inst.costs) for inst in target_instances[gt]]
    for gt in graph_types)

# Compute empirical homodist and P at target size (separate instances)
target_fit_instances = Dict(gt => generate_instances(gt, N_TARGET, NUM_FIT_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 50)) for gt in graph_types)
max_edges_target = maximum(inst.num_edges for gt in graph_types for inst in target_fit_instances[gt])

empirical_homodists_target = Dict{String, Array{Float64, 3}}()
sampled_homodists_target = Dict{String, Array{Float64, 3}}()
empirical_P_target = Dict{String, Vector{Float64}}()

println("Computing homodist at n=$N_TARGET (exact + sampled)...")
for gt in graph_types
    println("  $gt exact...")
    homodists_exact = map(target_fit_instances[gt]) do inst
        get_homogeneous_distribution_from_costs_direct(
            inst.costs, inst.num_edges, inst.num_vertices;
            max_num_edges=max_edges_target)
    end
    empirical_homodists_target[gt] = average_distributions(homodists_exact)

    println("  $gt sampled (S=$SAMPLES_PER_COST)...")
    homodists_sampled = map(target_fit_instances[gt]) do inst
        get_homogeneous_distribution_from_costs_sampled(
            inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
            max_num_edges=max_edges_target,
            rng=MersenneTwister(SEED + hash(gt) + hash(inst.num_edges)))
    end
    sampled_homodists_target[gt] = average_distributions(homodists_sampled)

    empirical_P_target[gt] = compute_empirical_P(target_fit_instances[gt], max_edges_target)
end

methods = ["Transfer", "PaperProxy", "Tri+TriP", "Norm+NormP", "EmpN+EmpP", "SampN+EmpP"]
results = Dict(gt => Dict(m => Dict{Int, Vector{Float64}}() for m in methods) for gt in graph_types)

for gt in graph_types
    println("\n--- $gt ---")

    for p in P_VALUES
        println("  p=$p:")

        # === TRANSFER ===
        if p == 1
            rng_opt = MersenneTwister(SEED + hash(gt) + p * 100)
            src_params = map(source_instances[gt]) do inst
                best_exp = -Inf; best_γ = 0.0; best_β = 0.0
                for γ in range(0.02, 2.0, length=80), β in range(0.01, π/2, length=80)
                    e = qaoa_expectation(inst.costs, N_SOURCE, [γ], [β])
                    if e > best_exp; best_exp = e; best_γ = γ; best_β = β; end
                end
                (γ=best_γ, β=best_β)
            end
            med_γ = [median([sp.γ for sp in src_params])]
            med_β = [median([sp.β for sp in src_params])]
            transfer_ratios = map(enumerate(target_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_TARGET, med_γ, med_β) / target_optimal[gt][i]
            end
        else
            src_ramp = map(source_instances[gt]) do inst
                best_params = (0.0, 0.0, 0.0, 0.0); best_exp = -Inf; gs = GRID_SIZE_RAMP
                for γ₁ in range(0.02, 0.40, length=gs),
                    γ_f in range(0.10, 0.70, length=gs),
                    β₁ in range(0.05, 0.45, length=gs),
                    β_f in range(0.01, 0.25, length=gs)
                    γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, p)
                    exp_val = qaoa_expectation(inst.costs, N_SOURCE, γs_pi .* π, βs_pi .* π)
                    if exp_val > best_exp; best_exp = exp_val; best_params = (γ₁, γ_f, β₁, β_f); end
                end
                best_params
            end
            med_ramp = Tuple(median([sp[i] for sp in src_ramp]) for i in 1:4)
            γs_t, βs_t = linear_ramp(med_ramp..., p)
            transfer_ratios = map(enumerate(target_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_TARGET, γs_t .* π, βs_t .* π) / target_optimal[gt][i]
            end
        end
        results[gt]["Transfer"][p] = transfer_ratios
        print_stats("Transfer", transfer_ratios)

        # === PAPERPROXY ===
        proxy_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, N_TARGET)
            proxy = PaperProxy(inst.num_edges, N_TARGET, p_eff)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_TARGET, p)
            qaoa_expectation(inst.costs, N_TARGET, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["PaperProxy"][p] = proxy_ratios
        print_stats("PaperProxy", proxy_ratios)

        # === FITTED TRIANGLE + TRIANGLE P (consistent) ===
        tri_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            fp = fitted_tri_params[gt]
            proxy = IntuitiveTriangleProxy(inst.num_edges, N_TARGET, fp...)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_TARGET, p)
            qaoa_expectation(inst.costs, N_TARGET, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["Tri+TriP"][p] = tri_ratios
        print_stats("Tri+TriP", tri_ratios)

        # === FITTED NORMAL + NORMAL P (consistent) ===
        norm_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            fp = fitted_norm_params[gt]
            proxy = NormalProxy(inst.num_edges, N_TARGET, fp...)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, N_TARGET, p)
            qaoa_expectation(inst.costs, N_TARGET, best_γ, best_β) / target_optimal[gt][i]
        end
        results[gt]["Norm+NormP"][p] = norm_ratios
        print_stats("Norm+NormP", norm_ratios)

        # === EMPIRICAL N + EMPIRICAL P ===
        emp_homodist = empirical_homodists_target[gt]
        P_full = empirical_P_target[gt]
        m_hd = size(emp_homodist, 1) - 1
        P_vals = P_full[1:min(m_hd+1, length(P_full))]
        if length(P_vals) < m_hd + 1
            P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
        end
        best_γ_emp, best_β_emp, _ = optimize_via_proxy(emp_homodist, P_vals, N_TARGET, p)
        emp_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            qaoa_expectation(inst.costs, N_TARGET, best_γ_emp, best_β_emp) / target_optimal[gt][i]
        end
        results[gt]["EmpN+EmpP"][p] = emp_ratios
        print_stats("EmpN+EmpP", emp_ratios)

        # === SAMPLED N + EMPIRICAL P ===
        samp_homodist = sampled_homodists_target[gt]
        m_hd_s = size(samp_homodist, 1) - 1
        P_vals_s = P_full[1:min(m_hd_s+1, length(P_full))]
        if length(P_vals_s) < m_hd_s + 1
            P_vals_s = vcat(P_vals_s, zeros(m_hd_s + 1 - length(P_vals_s)))
        end
        best_γ_samp, best_β_samp, _ = optimize_via_proxy(samp_homodist, P_vals_s, N_TARGET, p)
        samp_ratios = map(enumerate(target_instances[gt])) do (i, inst)
            qaoa_expectation(inst.costs, N_TARGET, best_γ_samp, best_β_samp) / target_optimal[gt][i]
        end
        results[gt]["SampN+EmpP"][p] = samp_ratios
        print_stats("SampN+EmpP", samp_ratios)
    end
end


#==============================================================================#
#   PLOTS                                                                      #
#==============================================================================#

println("\n--- Plotting ---")

method_colors = Dict(
    "Transfer" => :steelblue,
    "PaperProxy" => :coral,
    "Tri+TriP" => :mediumseagreen,
    "Norm+NormP" => :purple,
    "EmpN+EmpP" => :gold,
    "SampN+EmpP" => :mediumpurple,
)

fig = Figure(size=(500 * length(graph_types), 400 * length(P_VALUES)))
for (pi, p) in enumerate(P_VALUES)
    for (gi, gt) in enumerate(graph_types)
        ax = Axis(fig[pi, gi],
            xlabel="Method", ylabel="Approx Ratio",
            title="$gt  p=$p",
            xticklabelrotation=π/4)

        for (mi, method) in enumerate(methods)
            if !haskey(results[gt][method], p); continue; end
            vals = results[gt][method][p]
            boxplot!(ax, fill(Float64(mi), length(vals)), vals,
                color=method_colors[method], width=0.6)
        end
        ax.xticks = (collect(1:length(methods)), collect(methods))
    end
end

Label(fig[0, :],
    "Systematic Proxy Evaluation: n=$N_SOURCE→$N_TARGET, $NUM_TARGET_INSTANCES instances\nFitted at n=$N_FIT on $NUM_FIT_INSTANCES instances | Consistent N/P pairing",
    fontsize=14, font=:bold)
save_figure(fig, "systematic_proxy_evaluation.png")


#==============================================================================#
#   SUMMARY TABLE                                                              #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY: Mean Approximation Ratios")
println("=" ^ 80)

for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-6s", "Graph")
    for m in methods; @printf("  %-14s", m); end
    println()
    println("  " * "-" ^ (6 + 16 * length(methods)))
    for gt in graph_types
        @printf("  %-6s", gt)
        for m in methods
            if haskey(results[gt][m], p)
                @printf("  %.4f±%.4f", mean(results[gt][m][p]), std(results[gt][m][p]))
            else
                @printf("  %14s", "—")
            end
        end
        println()
    end
end

# Key comparison: Which proxy is best for non-ER?
println("\n\nKey: Best proxy method for non-ER graphs (excluding Transfer)")
for p in P_VALUES
    println("  p=$p:")
    for gt in ["BA", "WS"]
        proxy_methods = ["PaperProxy", "Tri+TriP", "Norm+NormP", "EmpN+EmpP", "SampN+EmpP"]
        means = [(m, mean(results[gt][m][p])) for m in proxy_methods if haskey(results[gt][m], p)]
        best = sort(means, by=x -> -x[2])
        @printf("    %s: %s (%.4f)", gt, best[1][1], best[1][2])
        @printf(" > %s (%.4f)", best[2][1], best[2][2])
        @printf(" > %s (%.4f)\n", best[3][1], best[3][2])
    end
end

println("\nDone!")
