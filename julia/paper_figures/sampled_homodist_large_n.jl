#=
sampled_homodist_large_n.jl — Test sampled homodist at n=14-20 for QAOA
parameter optimization, where exact homodist becomes impractical.

The key question: can sampled homodist + proxy optimization produce useful
QAOA parameters at scales where brute-force homodist is too expensive?

Compare against:
  - Transfer baseline (doesn't need homodist)
  - PaperProxy (analytical, no homodist needed)

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf
using Distributions: Binomial, pdf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VALUES = [14, 16, 18]       # target graph sizes
const N_SOURCE = 9                  # source for transfer
const NUM_HOMODIST_INSTANCES = 20   # instances for homodist averaging
const NUM_EVAL_INSTANCES = 10       # instances for QAOA evaluation
const SAMPLES_PER_COST = 20        # higher S for better accuracy at large n

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

const P_VALUES = [1, 3]
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const SEED = 42

#==============================================================================#
#                    HELPERS                                                    #
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

function print_stats(label, values)
    @printf("  %-50s  mean=%.4f  std=%.4f\n", label, mean(values), std(values))
end


#==============================================================================#
#   MAIN EXPERIMENT                                                            #
#==============================================================================#

println("=" ^ 80)
println("Sampled Homodist at Large n: QAOA Parameter Optimization")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]
methods = ["Transfer", "PaperProxy", "SampN+EmpP"]

all_results = Dict{String, Dict{String, Dict{Int, Dict{Int, Vector{Float64}}}}}()
for gt in graph_types
    all_results[gt] = Dict(m => Dict{Int, Dict{Int, Vector{Float64}}}() for m in methods)
end

# Generate source instances for transfer (once)
source_instances = Dict(gt => generate_instances(gt, N_SOURCE, NUM_EVAL_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 888)) for gt in graph_types)

for n_target in N_VALUES
    println("\n" * "=" ^ 40)
    println("  n = $n_target")
    println("=" ^ 40)

    # Generate eval and homodist instances
    eval_instances = Dict(gt => generate_instances(gt, n_target, NUM_EVAL_INSTANCES;
        rng=MersenneTwister(SEED + hash(gt) + n_target * 100)) for gt in graph_types)
    eval_optimal = Dict(gt => [maxcut_optimal(inst.costs) for inst in eval_instances[gt]]
        for gt in graph_types)
    homodist_pool = Dict(gt => generate_instances(gt, n_target, NUM_HOMODIST_INSTANCES;
        rng=MersenneTwister(SEED + hash(gt) + n_target * 200)) for gt in graph_types)

    for gt in graph_types
        println("\n--- $gt (n=$n_target) ---")

        for p in P_VALUES
            println("  p=$p:")

            # === TRANSFER ===
            if p == 1
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
                transfer_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                    qaoa_expectation(inst.costs, n_target, med_γ, med_β) / eval_optimal[gt][i]
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
                transfer_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                    qaoa_expectation(inst.costs, n_target, γs_t .* π, βs_t .* π) / eval_optimal[gt][i]
                end
            end
            if !haskey(all_results[gt]["Transfer"], n_target)
                all_results[gt]["Transfer"][n_target] = Dict{Int, Vector{Float64}}()
            end
            all_results[gt]["Transfer"][n_target][p] = transfer_ratios
            print_stats("Transfer", transfer_ratios)

            # === PAPERPROXY ===
            pp_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, n_target)
                proxy = PaperProxy(inst.num_edges, n_target, p_eff)
                homodist = cpu_compute_homodist(proxy)
                P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                best_γ, best_β, _ = optimize_via_proxy(homodist, P_vals, n_target, p)
                qaoa_expectation(inst.costs, n_target, best_γ, best_β) / eval_optimal[gt][i]
            end
            if !haskey(all_results[gt]["PaperProxy"], n_target)
                all_results[gt]["PaperProxy"][n_target] = Dict{Int, Vector{Float64}}()
            end
            all_results[gt]["PaperProxy"][n_target][p] = pp_ratios
            print_stats("PaperProxy", pp_ratios)

            # === SAMPLED HOMODIST ===
            println("    Computing sampled homodist (S=$SAMPLES_PER_COST)...")
            max_edges = maximum(inst.num_edges for inst in homodist_pool[gt])
            t_sample = @elapsed begin
                sampled_homodists = map(homodist_pool[gt]) do inst
                    get_homogeneous_distribution_from_costs_sampled(
                        inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
                        max_num_edges=max_edges,
                        rng=MersenneTwister(SEED + hash(gt) + hash(inst.num_edges) + n_target))
                end
                sampled_avg = average_distributions(sampled_homodists)
            end
            P_emp = compute_empirical_P(homodist_pool[gt], max_edges)
            m_hd = size(sampled_avg, 1) - 1
            P_vals = P_emp[1:min(m_hd+1, length(P_emp))]
            if length(P_vals) < m_hd + 1
                P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
            end

            best_γ_s, best_β_s, _ = optimize_via_proxy(sampled_avg, P_vals, n_target, p)
            samp_ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, n_target, best_γ_s, best_β_s) / eval_optimal[gt][i]
            end
            if !haskey(all_results[gt]["SampN+EmpP"], n_target)
                all_results[gt]["SampN+EmpP"][n_target] = Dict{Int, Vector{Float64}}()
            end
            all_results[gt]["SampN+EmpP"][n_target][p] = samp_ratios
            print_stats("SampN+EmpP ($(round(t_sample, digits=2))s)", samp_ratios)
        end
    end
end


#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY: Approximation Ratios at Large n")
println("=" ^ 80)

for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-5s %-4s", "Graph", "n")
    for m in methods; @printf("  %-16s", m); end
    println()
    println("  " * "-" ^ (10 + 18 * length(methods)))
    for gt in graph_types
        for n in N_VALUES
            @printf("  %-5s %-4d", gt, n)
            for m in methods
                if haskey(all_results[gt][m], n) && haskey(all_results[gt][m][n], p)
                    vals = all_results[gt][m][n][p]
                    @printf("  %.4f±%.4f  ", mean(vals), std(vals))
                else
                    @printf("  %-16s", "—")
                end
            end
            println()
        end
    end
end

println("\nDone!")
