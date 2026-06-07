#=
smoothed_homodist_investigation.jl — Test smoothed (regularized) homodist
for QAOA parameter optimization.

Instead of fitting a parametric proxy, this applies Gaussian smoothing to the
sampled empirical homodist. The hypothesis: smoothing removes high-frequency
noise that causes proxy landscape overfitting at high depth, while preserving
the important structural features.

Compare smoothed homodist at different smoothing levels against:
  - Raw SampN+EmpP (no smoothing)
  - Transfer baseline
  - PaperProxy

Also test: does PaperProxy succeed partly because it acts like a
maximally-smoothed version of the true homodist? If so, intermediate
smoothing levels should interpolate between SampN+EmpP and PaperProxy.

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VALUES = [12, 14, 16]
const N_SOURCE = 9
const NUM_HOMODIST_INSTANCES = 20
const NUM_EVAL_INSTANCES = 10
const SAMPLES_PER_COST = 20

const P_VALUES = [1, 3, 5]
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const SEED = 42

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

# Smoothing sigma values to test (in units of array indices)
const SMOOTH_SIGMAS = [0.5, 1.0, 2.0, 4.0]

#==============================================================================#
#                    SMOOTHING                                                  #
#==============================================================================#

"""
Apply 3D Gaussian smoothing to a homodist array.
Uses separable 1D convolutions along each axis for efficiency.
"""
function smooth_homodist(homodist::Array{Float64, 3}, sigma::Float64)
    if sigma <= 0
        return copy(homodist)
    end

    # Create 1D Gaussian kernel
    radius = ceil(Int, 3 * sigma)
    kernel = [exp(-0.5 * (x / sigma)^2) for x in -radius:radius]
    kernel ./= sum(kernel)

    result = copy(homodist)

    # Convolve along each dimension
    for dim in 1:3
        result = convolve_along_dim(result, kernel, dim)
    end

    # Ensure non-negative
    result .= max.(result, 0.0)

    return result
end

"""Convolve array along a specific dimension with a 1D kernel."""
function convolve_along_dim(arr::Array{Float64, 3}, kernel::Vector{Float64}, dim::Int)
    result = zeros(size(arr))
    radius = div(length(kernel), 2)
    sz = size(arr, dim)

    for I in CartesianIndices(arr)
        idx = Tuple(I)
        val = 0.0
        for k in -radius:radius
            j = idx[dim] + k
            # Reflect at boundaries
            j = clamp(j, 1, sz)
            new_idx = ntuple(d -> d == dim ? j : idx[d], 3)
            val += arr[new_idx...] * kernel[k + radius + 1]
        end
        result[I] = val
    end

    return result
end

"""
Smooth the P distribution (1D Gaussian smoothing).
"""
function smooth_P(P::Vector{Float64}, sigma::Float64)
    if sigma <= 0
        return copy(P)
    end
    radius = ceil(Int, 3 * sigma)
    kernel = [exp(-0.5 * (x / sigma)^2) for x in -radius:radius]
    kernel ./= sum(kernel)

    n = length(P)
    result = zeros(n)
    for i in 1:n
        for k in -radius:radius
            j = clamp(i + k, 1, n)
            result[i] += P[j] * kernel[k + radius + 1]
        end
    end
    return result
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
println("Smoothed Homodist Investigation")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]
base_methods = ["Random", "Transfer", "PaperProxy", "SampN+EmpP"]
smooth_methods = ["Smooth(σ=$(σ))" for σ in SMOOTH_SIGMAS]
all_method_names = vcat(base_methods, smooth_methods)

# results[gt][method][n][p] = Vector{Float64}
all_results = Dict{String, Dict{String, Dict{Int, Dict{Int, Vector{Float64}}}}}()
for gt in graph_types
    all_results[gt] = Dict(m => Dict{Int, Dict{Int, Vector{Float64}}}() for m in all_method_names)
end

for n_target in N_VALUES
    for gt in graph_types
        println("\n" * "=" ^ 40)
        println("  $gt (n=$n_target)")
        println("=" ^ 40)

        # Generate instances
        eval_rng = MersenneTwister(SEED + hash(gt) + n_target * 100)
        eval_insts = generate_instances(gt, n_target, NUM_EVAL_INSTANCES; rng=eval_rng)
        eval_optimal = [maxcut_optimal(inst.costs) for inst in eval_insts]

        hd_rng = MersenneTwister(SEED + hash(gt) + n_target * 200)
        hd_insts = generate_instances(gt, n_target, NUM_HOMODIST_INSTANCES; rng=hd_rng)

        src_rng = MersenneTwister(SEED + hash(gt) + 888)
        source_insts = generate_instances(gt, N_SOURCE, NUM_EVAL_INSTANCES; rng=src_rng)

        # Compute sampled homodist
        max_edges = maximum(inst.num_edges for inst in hd_insts)
        println("  Computing sampled homodist...")
        sampled_homodists = map(hd_insts) do inst
            get_homogeneous_distribution_from_costs_sampled(
                inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
                max_num_edges=max_edges,
                rng=MersenneTwister(SEED + hash(gt) + hash(inst.num_edges) + n_target))
        end
        raw_homodist = average_distributions(sampled_homodists)

        # Empirical P
        P_emp = compute_empirical_P(hd_insts, max_edges)
        m_hd = size(raw_homodist, 1) - 1
        P_vals = P_emp[1:min(m_hd+1, length(P_emp))]
        if length(P_vals) < m_hd + 1
            P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
        end

        # Pre-compute smoothed homodists and P values
        smoothed_hds = Dict{Float64, Array{Float64, 3}}()
        smoothed_Ps = Dict{Float64, Vector{Float64}}()
        for σ in SMOOTH_SIGMAS
            println("  Smoothing with σ=$σ...")
            smoothed_hds[σ] = smooth_homodist(raw_homodist, σ)
            smoothed_Ps[σ] = smooth_P(P_vals, σ)
        end

        for p in P_VALUES
            println("\n  p=$p:")

            # === RANDOM ===
            random_ratios = [mean(inst.costs) / eval_optimal[i] for (i, inst) in enumerate(eval_insts)]
            if !haskey(all_results[gt]["Random"], n_target)
                all_results[gt]["Random"][n_target] = Dict{Int, Vector{Float64}}()
            end
            all_results[gt]["Random"][n_target][p] = random_ratios
            print_stats("Random baseline", random_ratios)

            # === TRANSFER ===
            t_γ, t_β = transfer_optimize(source_insts, N_SOURCE, p)
            transfer_ratios = [qaoa_expectation(inst.costs, n_target, t_γ, t_β) / eval_optimal[i]
                               for (i, inst) in enumerate(eval_insts)]
            if !haskey(all_results[gt]["Transfer"], n_target)
                all_results[gt]["Transfer"][n_target] = Dict{Int, Vector{Float64}}()
            end
            all_results[gt]["Transfer"][n_target][p] = transfer_ratios
            print_stats("Transfer", transfer_ratios)

            # === PAPERPROXY ===
            pp_ratios = map(enumerate(eval_insts)) do (i, inst)
                p_eff = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, n_target)
                proxy = PaperProxy(inst.num_edges, n_target, p_eff)
                hd = cpu_compute_homodist(proxy)
                Pv = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                best_γ, best_β, _ = optimize_via_proxy(hd, Pv, n_target, p)
                qaoa_expectation(inst.costs, n_target, best_γ, best_β) / eval_optimal[i]
            end
            if !haskey(all_results[gt]["PaperProxy"], n_target)
                all_results[gt]["PaperProxy"][n_target] = Dict{Int, Vector{Float64}}()
            end
            all_results[gt]["PaperProxy"][n_target][p] = pp_ratios
            print_stats("PaperProxy", pp_ratios)

            # === RAW SAMPN+EMPP ===
            best_γ, best_β, _ = optimize_via_proxy(raw_homodist, P_vals, n_target, p)
            samp_ratios = [qaoa_expectation(inst.costs, n_target, best_γ, best_β) / eval_optimal[i]
                           for (i, inst) in enumerate(eval_insts)]
            if !haskey(all_results[gt]["SampN+EmpP"], n_target)
                all_results[gt]["SampN+EmpP"][n_target] = Dict{Int, Vector{Float64}}()
            end
            all_results[gt]["SampN+EmpP"][n_target][p] = samp_ratios
            print_stats("SampN+EmpP (raw)", samp_ratios)

            # === SMOOTHED VERSIONS ===
            for σ in SMOOTH_SIGMAS
                method_name = "Smooth(σ=$(σ))"
                best_γ_s, best_β_s, _ = optimize_via_proxy(smoothed_hds[σ], smoothed_Ps[σ], n_target, p)
                smooth_ratios = [qaoa_expectation(inst.costs, n_target, best_γ_s, best_β_s) / eval_optimal[i]
                                 for (i, inst) in enumerate(eval_insts)]
                if !haskey(all_results[gt][method_name], n_target)
                    all_results[gt][method_name][n_target] = Dict{Int, Vector{Float64}}()
                end
                all_results[gt][method_name][n_target][p] = smooth_ratios
                print_stats(method_name, smooth_ratios)
            end
        end
    end
end

#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-5s %-4s", "Graph", "n")
    for m in all_method_names; @printf("  %-14s", m); end
    println()
    println("  " * "-" ^ (12 + 16 * length(all_method_names)))
    for gt in graph_types
        for n in N_VALUES
            @printf("  %-5s %-4d", gt, n)
            for m in all_method_names
                if haskey(all_results[gt][m], n) && haskey(all_results[gt][m][n], p)
                    vals = all_results[gt][m][n][p]
                    @printf("  %.4f        ", mean(vals))
                else
                    @printf("  %-14s", "—")
                end
            end
            println()
        end
    end
end

# Special analysis: does smoothing help at high depth?
println("\n\nSMOOTHING EFFECT AT HIGH DEPTH:")
println("(Difference from raw SampN+EmpP, positive = smoothing helps)")
for p in P_VALUES
    println("\n  p=$p:")
    for gt in graph_types
        for n in N_VALUES
            raw_val = haskey(all_results[gt]["SampN+EmpP"], n) && haskey(all_results[gt]["SampN+EmpP"][n], p) ?
                mean(all_results[gt]["SampN+EmpP"][n][p]) : NaN
            if isnan(raw_val); continue; end
            @printf("  %-5s n=%-4d  raw=%.4f", gt, n, raw_val)
            for σ in SMOOTH_SIGMAS
                mn = "Smooth(σ=$(σ))"
                sv = haskey(all_results[gt][mn], n) && haskey(all_results[gt][mn][n], p) ?
                    mean(all_results[gt][mn][n][p]) : NaN
                if !isnan(sv)
                    @printf("  σ=%.1f: %+.4f", σ, sv - raw_val)
                end
            end
            println()
        end
    end
end

println("\nDone!")
