#=
computational_scaling.jl — Timing analysis for paper: how do different methods
scale with n?

Measures wall-clock time for each step:
  1. Cost computation (maxcut_costs): O(2^n * m)
  2. Sampled homodist: O(2^n * S * n * num_instances)
  3. Proxy optimization (grid search): O(m^2 * n * p * K_grid)
  4. Transfer source optimization: O(2^n_src * p * K_grid * num_src)
  5. QAOA evaluation (single): O(2^n * n * p)
  6. Coord descent refinement: O(2^n * n * p * steps)

This provides data for a paper figure showing computational feasibility.

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_VALUES = [10, 12, 14, 16, 18, 20, 22]
const N_SOURCE = 9
const NUM_INSTANCES = 10
const S = 20
const ER_P_EDGE = 0.5
const P_DEPTH = 3

const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const SEED = 42

#==============================================================================#
#                    HELPERS                                                    #
#==============================================================================#

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

#==============================================================================#
#   MAIN EXPERIMENT                                                            #
#==============================================================================#

println("=" ^ 80)
println("Computational Scaling Analysis (ER, p=$P_DEPTH)")
println("=" ^ 80)

# Warmup at small n
println("Warming up...")
warmup_inst = generate_er_instance(8, ER_P_EDGE; rng=MersenneTwister(0))
_ = qaoa_expectation(warmup_inst.costs, 8, [0.5], [0.3])
warmup_homodist = get_homogeneous_distribution_from_costs_sampled(
    warmup_inst.costs, warmup_inst.num_edges, 8, 5;
    rng=MersenneTwister(0))
proxy = PaperProxy(warmup_inst.num_edges, 8, 0.5)
_ = cpu_compute_homodist(proxy)

# Source transfer params (one-time cost)
source_instances = [generate_er_instance(N_SOURCE, ER_P_EDGE;
    rng=MersenneTwister(SEED + 888 + i)) for i in 1:10]

println("\nMeasuring transfer source optimization (one-time cost)...")
t_transfer_source = @elapsed begin
    gs = GRID_SIZE_RAMP
    src_ramp = map(source_instances) do inst
        best_params = (0.0, 0.0, 0.0, 0.0); best_exp = -Inf
        for γ₁ in range(0.02, 0.40, length=gs),
            γ_f in range(0.10, 0.70, length=gs),
            β₁ in range(0.05, 0.45, length=gs),
            β_f in range(0.01, 0.25, length=gs)
            γs_pi, βs_pi = linear_ramp(γ₁, γ_f, β₁, β_f, P_DEPTH)
            exp_val = qaoa_expectation(inst.costs, N_SOURCE, γs_pi .* π, βs_pi .* π)
            if exp_val > best_exp; best_exp = exp_val; best_params = (γ₁, γ_f, β₁, β_f); end
        end
        best_params
    end
    med_ramp = Tuple(median([sp[i] for sp in src_ramp]) for i in 1:4)
    transfer_γ, transfer_β = linear_ramp(med_ramp..., P_DEPTH)
    transfer_γ = transfer_γ .* π
    transfer_β = transfer_β .* π
end
@printf("Transfer source optimization: %.2fs (one-time, n_source=%d)\n", t_transfer_source, N_SOURCE)

# Per-n timing
timings = Dict{Int, Dict{String, Float64}}()

for n in N_VALUES
    println("\n--- n=$n (2^n=$(1 << n)) ---")
    timings[n] = Dict{String, Float64}()

    # Cost computation
    t_cost = @elapsed begin
        instances = [generate_er_instance(n, ER_P_EDGE;
            rng=MersenneTwister(SEED + n * 100 + i)) for i in 1:NUM_INSTANCES]
    end
    timings[n]["costs"] = t_cost
    @printf("  Cost computation (%d instances): %.3fs\n", NUM_INSTANCES, t_cost)

    max_edges = maximum(inst.num_edges for inst in instances)

    # Sampled homodist
    t_homodist = @elapsed begin
        sampled_homodists = map(instances) do inst
            get_homogeneous_distribution_from_costs_sampled(
                inst.costs, inst.num_edges, inst.num_vertices, S;
                max_num_edges=max_edges,
                rng=MersenneTwister(SEED + hash(inst.num_edges) + n))
        end
        sampled_avg = average_distributions(sampled_homodists)
    end
    timings[n]["homodist"] = t_homodist
    @printf("  Sampled homodist (S=%d, %d inst): %.3fs\n", S, NUM_INSTANCES, t_homodist)

    P_emp = compute_empirical_P(instances, max_edges)
    m_hd = size(sampled_avg, 1) - 1
    P_vals = P_emp[1:min(m_hd+1, length(P_emp))]
    if length(P_vals) < m_hd + 1
        P_vals = vcat(P_vals, zeros(m_hd + 1 - length(P_vals)))
    end

    # Proxy optimization
    t_proxy = @elapsed begin
        proxy_γ, proxy_β, _ = optimize_via_proxy(sampled_avg, P_vals, n, P_DEPTH)
    end
    timings[n]["proxy_opt"] = t_proxy
    @printf("  Proxy optimization (10^4 ramps): %.3fs\n", t_proxy)

    # PaperProxy homodist (analytical)
    t_pp = @elapsed begin
        pp = PaperProxy(instances[1].num_edges, n, ER_P_EDGE)
        pp_homodist = cpu_compute_homodist(pp)
        pp_P = [P_cost_distribution(pp, c) for c in 0:instances[1].num_edges]
    end
    timings[n]["paperproxy"] = t_pp
    @printf("  PaperProxy homodist: %.3fs\n", t_pp)

    # QAOA evaluation (single instance)
    t_eval = @elapsed begin
        _ = qaoa_expectation(instances[1].costs, n, proxy_γ, proxy_β)
    end
    timings[n]["qaoa_eval"] = t_eval
    @printf("  QAOA eval (single): %.4fs\n", t_eval)

    # Coord descent refinement (single instance)
    if n <= 20  # too slow at n=22
        t_refine = @elapsed begin
            _ = refine_qaoa_from(instances[1].costs, n, proxy_γ, proxy_β)
        end
        timings[n]["refine"] = t_refine
        @printf("  Coord descent refine (1 restart): %.2fs\n", t_refine)
    end

    # Total SampN+EmpP workflow
    total_samp = t_cost + t_homodist + t_proxy + t_eval
    timings[n]["total_sampn"] = total_samp
    @printf("  Total SampN+EmpP: %.2fs\n", total_samp)
end

#==============================================================================#
#   SUMMARY TABLE                                                              #
#==============================================================================#

println("\n" * "=" ^ 80)
println("TIMING SUMMARY (seconds, p=$P_DEPTH)")
println("=" ^ 80)

@printf("%-6s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n",
    "n", "Costs", "SampHD", "ProxyOpt", "PaperPx", "QAOAeval", "Refine", "TotalSamp")
println("-" ^ 76)
for n in N_VALUES
    t = timings[n]
    refine_str = haskey(t, "refine") ? @sprintf("%.3f", t["refine"]) : "N/A"
    @printf("%-6d %-10.3f %-10.3f %-10.3f %-10.3f %-10.4f %-10s %-10.2f\n",
        n, t["costs"], t["homodist"], t["proxy_opt"], t["paperproxy"],
        t["qaoa_eval"], refine_str, t["total_sampn"])
end

println("\nKey insight: SampN+EmpP total cost is dominated by:")
println("  - At small n: proxy optimization (grid search)")
println("  - At large n: cost computation + sampled homodist")
println("  Transfer cost is dominated by source optimization ($(round(t_transfer_source, digits=1))s one-time)")
println("  but requires no per-instance computation at deployment time.")

println("\nDone!")
