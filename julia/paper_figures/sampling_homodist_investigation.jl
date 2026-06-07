#=
sampling_homodist_investigation.jl — Validate sampling-based homodist estimation.

Questions answered:
  1. How does sampled N(c';d,c) compare to exact N(c';d,c)?
     (MSE, Pearson correlation vs sample count)
  2. How many samples per cost class are needed for the sampled homodist
     to produce proxy-optimal parameters that perform well on real QAOA?
  3. Does sampling work across graph types (ER, BA, WS)?
  4. How does computation time scale with sample count?

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf
using Distributions: Binomial, pdf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VERTICES = 12
const NUM_HOMODIST_INSTANCES = 30   # instances to average homodist over
const NUM_EVAL_INSTANCES = 20       # instances to evaluate QAOA on
const SAMPLES_PER_COST_VALUES = [1, 2, 5, 10, 20, 50, 100]
const NUM_SAMPLING_TRIALS = 5       # repeat sampling to measure variance

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

const P_VALUES = [1, 3]
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const SEED = 123

#==============================================================================#
#                    HELPERS                                                    #
#==============================================================================#

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

function print_stats(label, values)
    @printf("  %-45s  mean=%.4f  std=%.4f\n", label, mean(values), std(values))
end


#==============================================================================#
#   EXPERIMENT 1: Homodist Accuracy vs Sample Count                            #
#==============================================================================#

println("=" ^ 80)
println("EXPERIMENT 1: Sampled Homodist Accuracy vs Exact")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]

# Generate homodist instances
homodist_instances = Dict(gt => generate_instances(gt, N_VERTICES, NUM_HOMODIST_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt))) for gt in graph_types)

max_edges = maximum(inst.num_edges for gt in graph_types for inst in homodist_instances[gt])

# Compute exact homodist (ground truth)
println("\nComputing exact homodist (ground truth)...")
exact_homodists = Dict{String, Array{Float64, 3}}()
exact_P = Dict{String, Vector{Float64}}()

for gt in graph_types
    println("  $gt...")
    t = @elapsed begin
        homodists = map(homodist_instances[gt]) do inst
            get_homogeneous_distribution_from_costs_direct(
                inst.costs, inst.num_edges, inst.num_vertices;
                max_num_edges=max_edges)
        end
        exact_homodists[gt] = average_distributions(homodists)
    end
    exact_P[gt] = compute_empirical_P(homodist_instances[gt], max_edges)
    @printf("  %s: %.2fs\n", gt, t)
end

# Compute sampled homodist at various sample counts
println("\nComputing sampled homodist at various sample counts...")

# Store results: graph_type -> samples_per_cost -> trial -> (mse, pearson, time)
sampling_results = Dict{String, Dict{Int, Vector{NamedTuple{(:mse, :pearson_mean, :time), Tuple{Float64, Float64, Float64}}}}}()

for gt in graph_types
    println("\n--- $gt ---")
    sampling_results[gt] = Dict()

    for spc in SAMPLES_PER_COST_VALUES
        sampling_results[gt][spc] = []

        for trial in 1:NUM_SAMPLING_TRIALS
            trial_rng = MersenneTwister(SEED + hash(gt) + spc * 1000 + trial)
            t = @elapsed begin
                sampled_homodists = map(homodist_instances[gt]) do inst
                    get_homogeneous_distribution_from_costs_sampled(
                        inst.costs, inst.num_edges, inst.num_vertices, spc;
                        max_num_edges=max_edges, rng=MersenneTwister(rand(trial_rng, UInt64)))
                end
                sampled_avg = average_distributions(sampled_homodists)
            end

            # MSE between sampled and exact
            exact = exact_homodists[gt]
            sampled_avg_matched, exact_matched = pad_to_match(sampled_avg, exact)
            mse = sum((sampled_avg_matched .- exact_matched).^2) / length(exact_matched)

            # Pearson correlation per cost class
            pearsons = get_pearson_correlation_coefficients(sampled_avg_matched, exact_matched)
            # Filter out NaN (empty cost classes)
            valid_pearsons = filter(!isnan, pearsons)
            pearson_mean = isempty(valid_pearsons) ? NaN : mean(valid_pearsons)

            push!(sampling_results[gt][spc], (mse=mse, pearson_mean=pearson_mean, time=t))
        end

        mean_mse = mean(r.mse for r in sampling_results[gt][spc])
        mean_pearson = mean(r.pearson_mean for r in sampling_results[gt][spc])
        mean_time = mean(r.time for r in sampling_results[gt][spc])
        @printf("  S=%3d: MSE=%.2e  Pearson=%.6f  time=%.2fs\n", spc, mean_mse, mean_pearson, mean_time)
    end
end


#==============================================================================#
#   EXPERIMENT 2: QAOA Performance with Sampled Homodist                       #
#==============================================================================#

println("\n" * "=" ^ 80)
println("EXPERIMENT 2: QAOA Approximation Ratio with Sampled vs Exact Homodist")
println("=" ^ 80)

# Generate evaluation instances (separate from homodist instances)
eval_instances = Dict(gt => generate_instances(gt, N_VERTICES, NUM_EVAL_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 999)) for gt in graph_types)
eval_optimal = Dict(gt => [maxcut_optimal(inst.costs) for inst in eval_instances[gt]]
    for gt in graph_types)

# Also generate source instances for transfer baseline
const N_SMALL = 9
source_instances = Dict(gt => generate_instances(gt, N_SMALL, NUM_EVAL_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 888)) for gt in graph_types)

# Methods to test
test_samples = [2, 5, 10, 20, 50]
methods = vcat(["Transfer", "Exact"], ["S=$s" for s in test_samples])

qaoa_results = Dict(gt => Dict(m => Dict{Int, Vector{Float64}}() for m in methods) for gt in graph_types)

for gt in graph_types
    println("\n--- $gt ---")

    for p in P_VALUES
        println("  p=$p:")

        # === TRANSFER BASELINE ===
        if p == 1
            # Optimize on source, transfer median to target
            rng_opt = MersenneTwister(SEED + hash(gt) + p * 100)
            src_params = map(source_instances[gt]) do inst
                best_exp = -Inf; best_γ = 0.0; best_β = 0.0
                for γ in range(0.02, 2.0, length=80), β in range(0.01, π/2, length=80)
                    e = qaoa_expectation(inst.costs, N_SMALL, [γ], [β])
                    if e > best_exp; best_exp = e; best_γ = γ; best_β = β; end
                end
                (γ=best_γ, β=best_β)
            end
            med_γ = [median([sp.γ for sp in src_params])]
            med_β = [median([sp.β for sp in src_params])]
            transfer_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_VERTICES, med_γ, med_β) / eval_optimal[gt][i]
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
            transfer_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_VERTICES, γs_t .* π, βs_t .* π) / eval_optimal[gt][i]
            end
        end
        qaoa_results[gt]["Transfer"][p] = transfer_ratios
        print_stats("Transfer", transfer_ratios)

        # === EXACT HOMODIST (EmpN+EmpP) ===
        homodist_exact = exact_homodists[gt]
        P_exact = exact_P[gt]
        m_hd = size(homodist_exact, 1) - 1
        P_vals = P_exact[1:min(m_hd+1, length(P_exact))]
        if length(P_vals) < m_hd + 1
            P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
        end

        best_γ_exact, best_β_exact, _ = optimize_via_proxy(homodist_exact, P_vals, N_VERTICES, p)
        exact_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
            qaoa_expectation(inst.costs, N_VERTICES, best_γ_exact, best_β_exact) / eval_optimal[gt][i]
        end
        qaoa_results[gt]["Exact"][p] = exact_ratios
        print_stats("Exact homodist", exact_ratios)

        # === SAMPLED HOMODIST at various sample counts ===
        for spc in test_samples
            trial_rng = MersenneTwister(SEED + hash(gt) + spc * 2000 + p)
            sampled_homodists = map(homodist_instances[gt]) do inst
                get_homogeneous_distribution_from_costs_sampled(
                    inst.costs, inst.num_edges, inst.num_vertices, spc;
                    max_num_edges=max_edges, rng=MersenneTwister(rand(trial_rng, UInt64)))
            end
            sampled_avg = average_distributions(sampled_homodists)

            # Compute empirical P from the same instances
            P_sampled = compute_empirical_P(homodist_instances[gt], max_edges)
            m_hd_s = size(sampled_avg, 1) - 1
            P_vals_s = P_sampled[1:min(m_hd_s+1, length(P_sampled))]
            if length(P_vals_s) < m_hd_s + 1
                P_vals_s = vcat(P_vals_s, zeros(m_hd_s + 1 - length(P_vals_s)))
            end

            best_γ_s, best_β_s, _ = optimize_via_proxy(sampled_avg, P_vals_s, N_VERTICES, p)
            sampled_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_VERTICES, best_γ_s, best_β_s) / eval_optimal[gt][i]
            end
            qaoa_results[gt]["S=$spc"][p] = sampled_ratios
            print_stats("Sampled (S=$spc)", sampled_ratios)
        end
    end
end


#==============================================================================#
#   EXPERIMENT 3: Timing Comparison                                            #
#==============================================================================#

println("\n" * "=" ^ 80)
println("EXPERIMENT 3: Computation Time (single instance)")
println("=" ^ 80)

for gt in graph_types
    inst = homodist_instances[gt][1]
    println("\n--- $gt (n=$N_VERTICES, m=$(inst.num_edges)) ---")

    # Exact
    t_exact = @elapsed get_homogeneous_distribution_from_costs_direct(
        inst.costs, inst.num_edges, inst.num_vertices)
    @printf("  Exact:      %.4fs\n", t_exact)

    # Sampled at various counts
    for spc in [1, 5, 10, 20, 50, 100]
        t_sampled = @elapsed get_homogeneous_distribution_from_costs_sampled(
            inst.costs, inst.num_edges, inst.num_vertices, spc;
            rng=MersenneTwister(42))
        speedup = t_exact / max(t_sampled, 1e-6)
        @printf("  S=%-3d:      %.4fs  (%.1fx speedup)\n", spc, t_sampled, speedup)
    end
end


#==============================================================================#
#   PLOTS                                                                      #
#==============================================================================#

println("\n--- Plotting ---")

# Plot 1: MSE vs sample count
fig1 = Figure(size=(900, 400))
ax1 = Axis(fig1[1, 1],
    xlabel="Samples per cost class", ylabel="MSE vs exact",
    title="Homodist Accuracy (n=$N_VERTICES, averaged over $NUM_HOMODIST_INSTANCES instances)",
    xscale=log10, yscale=log10)

gt_colors = Dict("ER" => :steelblue, "BA" => :coral, "WS" => :mediumseagreen)
for gt in graph_types
    xs = Float64.(SAMPLES_PER_COST_VALUES)
    ys = [mean(r.mse for r in sampling_results[gt][s]) for s in SAMPLES_PER_COST_VALUES]
    yerr_lo = [mean(r.mse for r in sampling_results[gt][s]) - minimum(r.mse for r in sampling_results[gt][s]) for s in SAMPLES_PER_COST_VALUES]
    yerr_hi = [maximum(r.mse for r in sampling_results[gt][s]) - mean(r.mse for r in sampling_results[gt][s]) for s in SAMPLES_PER_COST_VALUES]
    scatterlines!(ax1, xs, ys, label=gt, color=gt_colors[gt], markersize=8)
    errorbars!(ax1, xs, ys, yerr_lo, yerr_hi, color=gt_colors[gt], whiskerwidth=6)
end
axislegend(ax1)
save_figure(fig1, "sampling_homodist_mse.png")

# Plot 2: QAOA approximation ratio vs sample count
fig2 = Figure(size=(500 * length(graph_types), 400 * length(P_VALUES)))
for (pi, p) in enumerate(P_VALUES)
    for (gi, gt) in enumerate(graph_types)
        ax = Axis(fig2[pi, gi],
            xlabel="Method", ylabel="Approx Ratio",
            title="$gt, p=$p",
            xticklabelrotation=π/4)

        method_list = methods
        positions = 1:length(method_list)
        for (mi, method) in enumerate(method_list)
            if !haskey(qaoa_results[gt][method], p); continue; end
            vals = qaoa_results[gt][method][p]
            boxplot!(ax, fill(Float64(mi), length(vals)), vals,
                color=mi <= 2 ? (:steelblue, :coral)[mi] : :mediumpurple,
                width=0.6)
        end
        ax.xticks = (collect(positions), collect(method_list))
    end
end
Label(fig2[0, :],
    "QAOA Performance: Sampled Homodist vs Exact (n=$N_VERTICES)",
    fontsize=14, font=:bold)
save_figure(fig2, "sampling_homodist_qaoa.png")

# Plot 3: Pearson correlation vs sample count
fig3 = Figure(size=(900, 400))
ax3 = Axis(fig3[1, 1],
    xlabel="Samples per cost class", ylabel="Mean Pearson r vs exact",
    title="Homodist Correlation with Exact (n=$N_VERTICES)",
    xscale=log10)

for gt in graph_types
    xs = Float64.(SAMPLES_PER_COST_VALUES)
    ys = [mean(r.pearson_mean for r in sampling_results[gt][s]) for s in SAMPLES_PER_COST_VALUES]
    scatterlines!(ax3, xs, ys, label=gt, color=gt_colors[gt], markersize=8)
end
hlines!(ax3, [1.0], color=:gray, linestyle=:dash, label="Perfect")
axislegend(ax3, position=:rb)
save_figure(fig3, "sampling_homodist_pearson.png")


#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

println("\n1. Homodist Accuracy (MSE):")
@printf("  %-5s", "Graph")
for s in SAMPLES_PER_COST_VALUES; @printf("  S=%-4d", s); end
println()
for gt in graph_types
    @printf("  %-5s", gt)
    for s in SAMPLES_PER_COST_VALUES
        @printf("  %.1e", mean(r.mse for r in sampling_results[gt][s]))
    end
    println()
end

println("\n2. QAOA Approximation Ratios:")
for p in P_VALUES
    println("  p=$p:")
    @printf("  %-6s", "Graph")
    for m in methods; @printf("  %-10s", m); end
    println()
    for gt in graph_types
        @printf("  %-6s", gt)
        for m in methods
            if haskey(qaoa_results[gt][m], p)
                @printf("  %.4f±%.3f", mean(qaoa_results[gt][m][p]), std(qaoa_results[gt][m][p]))
            else
                @printf("  %10s", "—")
            end
        end
        println()
    end
end

println("\nDone!")
