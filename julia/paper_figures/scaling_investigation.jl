#=
scaling_investigation.jl — Scale non-ER proxy evaluation to larger graphs.

Tests Transfer vs PaperProxy(p_eff) vs EmpiricalN+EmpiricalP at n=14 and n=16.
Also profiles the cost of computing empirical N at each size.

Key question: Does the EmpN+EmpP approach scale? Computing N is O(2^(2n)),
so n=14 → ~268M ops, n=16 → ~4B ops. Need to verify computational feasibility
and whether the approach still works at larger n.

For transfer, source graphs are n-2 (n=12 or n=14).

Started: 2026-04-01
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

# Test sizes: (source_n, target_n, num_fit_instances)
const SIZE_CONFIGS = [
    (src=9, tgt=12, nfit=50),    # baseline (matches previous experiments)
    (src=12, tgt=14, nfit=30),   # medium
    (src=12, tgt=16, nfit=10),   # large (fewer instances due to cost)
]

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

const NUM_EVAL_INSTANCES = 20
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
#   MAIN EXPERIMENT                                                            #
#==============================================================================#

graph_types = ["ER", "BA", "WS"]
methods = ["Transfer", "PaperProxy", "EmpN+EmpP"]

# Results: results[size_label][gt][method][p] = Vector{Float64}
all_results = Dict{String, Dict{String, Dict{String, Dict{Int, Vector{Float64}}}}}()
timing_results = Dict{String, Dict{String, Float64}}()  # size -> gt -> seconds per instance

for cfg in SIZE_CONFIGS
    size_label = "n=$(cfg.src)→$(cfg.tgt)"
    println("\n" * "=" ^ 80)
    println("SIZE: $size_label ($(cfg.nfit) fitting instances)")
    println("=" ^ 80)

    all_results[size_label] = Dict()
    timing_results[size_label] = Dict{String, Float64}()

    for gt in graph_types
        println("\n--- $gt ---")
        all_results[size_label][gt] = Dict(m => Dict{Int, Vector{Float64}}() for m in methods)

        # Generate evaluation instances
        rng_src = MersenneTwister(SEED + hash(gt) + hash(cfg.tgt) + 1)
        rng_tgt = MersenneTwister(SEED + hash(gt) + hash(cfg.tgt) + 2)
        source = generate_instances(gt, cfg.src, NUM_EVAL_INSTANCES; rng=rng_src)
        target = generate_instances(gt, cfg.tgt, NUM_EVAL_INSTANCES; rng=rng_tgt)
        target_opt = [maxcut_optimal(inst.costs) for inst in target]

        avg_tgt_edges = mean(inst.num_edges for inst in target)
        @printf("  tgt_edges: %.1f  c*: %.1f±%.1f\n",
            avg_tgt_edges, mean(target_opt), std(target_opt))

        # Generate fitting instances and compute empirical homodist
        rng_fit = MersenneTwister(SEED + hash(gt) + hash(cfg.tgt) + 50)
        fit_insts = generate_instances(gt, cfg.tgt, cfg.nfit; rng=rng_fit)
        max_edges = maximum(inst.num_edges for inst in fit_insts)

        println("  Computing empirical homodist ($(cfg.nfit) instances at n=$(cfg.tgt))...")
        t_start = time()
        homodists = map(fit_insts) do inst
            get_homogeneous_distribution_from_costs_direct(
                inst.costs, inst.num_edges, inst.num_vertices;
                max_num_edges=max_edges)
        end
        t_elapsed = time() - t_start
        timing_results[size_label][gt] = t_elapsed / cfg.nfit
        @printf("  Homodist computation: %.1fs total (%.2fs/instance)\n", t_elapsed, t_elapsed / cfg.nfit)

        emp_homodist = average_distributions(homodists)
        emp_P = compute_empirical_P(fit_insts, max_edges)

        for p in P_VALUES
            println("  p=$p:")

            # === TRANSFER ===
            if p == 1
                rng_opt = MersenneTwister(SEED + hash(gt) + hash(cfg.tgt) + p * 100)
                src_params = map(source) do inst
                    γs, βs, _ = optimize_qaoa_cd(inst.costs, cfg.src, 1, N_RESTARTS; rng=rng_opt)
                    (γs=γs, βs=βs)
                end
                med_γ = [median([sp.γs[1] for sp in src_params])]
                med_β = [median([sp.βs[1] for sp in src_params])]
                transfer_ratios = map(enumerate(target)) do (i, inst)
                    qaoa_expectation(inst.costs, cfg.tgt, med_γ, med_β) / target_opt[i]
                end
            else
                src_ramp = map(source) do inst
                    best_params = (0.0, 0.0, 0.0, 0.0); best_exp = -Inf
                    gs = GRID_SIZE_RAMP
                    for γ₁ in range(0.02, 0.40, length=gs),
                        γ_f in range(0.10, 0.70, length=gs),
                        β₁ in range(0.05, 0.45, length=gs),
                        β_f in range(0.01, 0.25, length=gs)
                        γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, p)
                        exp_val = qaoa_expectation(inst.costs, cfg.src, γs_pi .* π, βs_pi .* π)
                        if exp_val > best_exp; best_exp = exp_val; best_params = (γ₁, γ_f, β₁, β_f); end
                    end
                    best_params
                end
                med_ramp = Tuple(median([sp[i] for sp in src_ramp]) for i in 1:4)
                γs_t, βs_t = linear_ramp(med_ramp..., p)
                transfer_ratios = map(enumerate(target)) do (i, inst)
                    qaoa_expectation(inst.costs, cfg.tgt, γs_t .* π, βs_t .* π) / target_opt[i]
                end
            end
            all_results[size_label][gt]["Transfer"][p] = transfer_ratios
            print_stats("Transfer", transfer_ratios)

            # === PAPERPROXY ===
            proxy_ratios = map(enumerate(target)) do (i, inst)
                p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, cfg.tgt)
                proxy = PaperProxy(inst.num_edges, cfg.tgt, p_eff)
                homodist = cpu_compute_homodist(proxy)
                P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, cfg.tgt, p)
                qaoa_expectation(inst.costs, cfg.tgt, best_γ, best_β) / target_opt[i]
            end
            all_results[size_label][gt]["PaperProxy"][p] = proxy_ratios
            print_stats("PaperProxy", proxy_ratios)

            # === EMPIRICAL N + EMPIRICAL P ===
            emp_ratios = map(enumerate(target)) do (i, inst)
                m_hd = size(emp_homodist, 1) - 1
                P_vals = emp_P[1:min(m_hd+1, length(emp_P))]
                if length(P_vals) < m_hd + 1
                    P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
                end
                best_γ, best_β, _ = optimize_via_proxy(emp_homodist, P_vals, cfg.tgt, p)
                qaoa_expectation(inst.costs, cfg.tgt, best_γ, best_β) / target_opt[i]
            end
            all_results[size_label][gt]["EmpN+EmpP"][p] = emp_ratios
            print_stats("EmpN+EmpP", emp_ratios)
        end
    end
end


#==============================================================================#
#   PLOTS                                                                      #
#==============================================================================#

println("\n--- Plotting ---")

method_colors = Dict(
    "Transfer" => :steelblue,
    "PaperProxy" => :coral,
    "EmpN+EmpP" => :gold)

# One figure per graph type, showing scaling
for gt in graph_types
    size_labels = [cfg for cfg in ["n=9→12", "n=12→14", "n=12→16"]
                   if haskey(all_results, cfg) && haskey(all_results[cfg], gt)]

    fig = Figure(size=(400 * length(size_labels), 500))
    for (si, sl) in enumerate(size_labels)
        ax = Axis(fig[1, si],
            xlabel="QAOA Depth p", ylabel="Approx Ratio",
            title="$gt $sl", xticks=P_VALUES)
        for p in P_VALUES
            nm = length(methods); tw = 0.7; sw = tw / nm
            for (mi, method) in enumerate(methods)
                if !haskey(all_results[sl][gt][method], p); continue; end
                offset = (mi - (nm + 1) / 2) * sw
                vals = all_results[sl][gt][method][p]
                boxplot!(ax, fill(Float64(p) + offset, length(vals)), vals,
                    color=method_colors[method], width=sw * 0.8,
                    label=(p == 1 && si == 1 ? method : nothing))
            end
        end
        if si == 1; axislegend(ax, position=:rb); end
    end
    Label(fig[0, :], "$gt: Scaling Comparison ($NUM_EVAL_INSTANCES instances)",
        fontsize=14, font=:bold)
    save_figure(fig, "scaling_$(lowercase(gt)).png")
end


#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

println("\nHomodist Computation Time (seconds per instance):")
@printf("  %-15s", "Size")
for gt in graph_types; @printf("  %-8s", gt); end
println()
println("  " * "-" ^ (15 + 10 * length(graph_types)))
for cfg in SIZE_CONFIGS
    sl = "n=$(cfg.src)→$(cfg.tgt)"
    @printf("  %-15s", sl)
    for gt in graph_types
        if haskey(timing_results[sl], gt)
            @printf("  %-8.2f", timing_results[sl][gt])
        else
            @printf("  %-8s", "N/A")
        end
    end
    println()
end

for sl in sort(collect(keys(all_results)))
    println("\n$sl:")
    for p in P_VALUES
        println("  p=$p:")
        @printf("  %-8s", "Graph")
        for m in methods; @printf("  %-16s", m); end
        println()
        println("  " * "-" ^ (8 + 18 * length(methods)))
        for gt in graph_types
            @printf("  %-8s", gt)
            for m in methods
                if haskey(all_results[sl][gt][m], p)
                    vals = all_results[sl][gt][m][p]
                    @printf("  %.4f±%.4f  ", mean(vals), std(vals))
                else
                    @printf("  %-16s", "N/A")
                end
            end
            println()
        end
    end
end

# Gap analysis: Transfer - EmpN+EmpP across sizes
println("\n\nTransfer - EmpN+EmpP Gap (positive = Transfer better):")
for sl in sort(collect(keys(all_results)))
    println("  $sl:")
    for p in P_VALUES
        @printf("    p=%d:", p)
        for gt in ["BA", "WS"]
            t = all_results[sl][gt]["Transfer"][p]
            e = all_results[sl][gt]["EmpN+EmpP"][p]
            @printf("  %s=%+.4f", gt, mean(t) - mean(e))
        end
        println()
    end
end

println("\nDone!")
