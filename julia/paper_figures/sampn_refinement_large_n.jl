#=
sampn_refinement_large_n.jl — Test SampN+EmpP + Refinement vs Transfer + Refinement
at large n on ER and non-ER graphs.

The hybrid warmstart investigation (n=12) showed Proxy+Refine beats Transfer by
0.01-0.03. But that used PaperProxy. SampN+EmpP is better than PaperProxy on ER
at large n. Does SampN+EmpP + Refine beat Transfer + Refine?

Methods:
  1. Transfer (no refine) — baseline
  2. Transfer + Refine — transfer params + coord descent on target
  3. SampN+EmpP (no refine) — proxy only
  4. SampN+EmpP + Refine — proxy warmstart + coord descent on target
  5. Random + Refine — control

At n=14-16, QAOA eval is O(2^n) per call. Coord descent does ~5 scales × 2p params
× 6 deltas = ~60p evals per restart. With 3 restarts, that's ~180p evals.
At n=16, p=3: 180*3 * 2^16 ≈ 35M ops — very manageable.

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_SOURCE = 9
const N_TARGETS = [14, 16, 18]
const NUM_HOMODIST_INSTANCES = 10
const NUM_EVAL_INSTANCES = 10
const SAMPLES_PER_COST = 20

const ER_P_EDGE = 0.5
const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3

const P_VALUES = [1, 3, 5]
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10

# Refinement parameters
const REFINE_RESTARTS = 3
const REFINE_PERTURBATION = 0.15
const RANDOM_RESTARTS = 10

const SEED = 42

#==============================================================================#
#                    HELPERS                                                    #
#==============================================================================#

"""Coordinate descent refinement from a specific starting point."""
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

"""Warmstart + perturbed restarts, refined on real QAOA."""
function hybrid_optimize(costs, n, init_γs, init_βs; n_restarts=REFINE_RESTARTS,
    perturb=REFINE_PERTURBATION, rng=Random.default_rng())
    p = length(init_γs)
    best_γs, best_βs, best_exp = refine_qaoa_from(costs, n, init_γs, init_βs)

    for _ in 1:n_restarts
        trial_γs = init_γs .+ randn(rng, p) .* perturb
        trial_βs = init_βs .+ randn(rng, p) .* perturb
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

function get_transfer_params(source_instances, p)
    n_src = source_instances[1].num_vertices
    if p == 1
        src_params = map(source_instances) do inst
            best_exp = -Inf; best_γ = 0.0; best_β = 0.0
            for γ in range(0.02, 2.0, length=80), β in range(0.01, π/2, length=80)
                e = qaoa_expectation(inst.costs, n_src, [γ], [β])
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
                exp_val = qaoa_expectation(inst.costs, n_src, γs_pi .* π, βs_pi .* π)
                if exp_val > best_exp; best_exp = exp_val; best_params = (γ₁, γ_f, β₁, β_f); end
            end
            best_params
        end
        med_ramp = Tuple(median([sp[i] for sp in src_ramp]) for i in 1:4)
        γs_t, βs_t = linear_ramp(med_ramp..., p)
        return γs_t .* π, βs_t .* π
    end
end

#==============================================================================#
#   MAIN EXPERIMENT                                                            #
#==============================================================================#

println("=" ^ 80)
println("SampN+EmpP + Refinement vs Transfer + Refinement")
println("Source n=$N_SOURCE, Targets n=$N_TARGETS")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]
methods = ["Transfer", "Transfer+Refine", "SampN+EmpP", "SampN+Refine", "Random+Refine"]

all_results = Dict{Tuple{String,Int,Int}, Dict{String, Vector{Float64}}}()

for gt in graph_types
    println("\n" * "#" ^ 80)
    println("  Graph type: $gt")
    println("#" ^ 80)

    # Source instances
    source_insts = generate_instances(gt, N_SOURCE, 10; rng=MersenneTwister(SEED + hash(gt) + 1))

    for n_target in N_TARGETS
        println("\n  n=$n_target:")

        # Target instances
        eval_insts = generate_instances(gt, n_target, NUM_EVAL_INSTANCES;
            rng=MersenneTwister(SEED + hash(gt) + n_target * 100))
        eval_optimal = [maxcut_optimal(inst.costs) for inst in eval_insts]

        # Homodist pool (for SampN+EmpP)
        homodist_pool = generate_instances(gt, n_target, NUM_HOMODIST_INSTANCES;
            rng=MersenneTwister(SEED + hash(gt) + n_target * 200))

        # Precompute sampled homodist
        max_edges = maximum(inst.num_edges for inst in homodist_pool)
        sampled_homodists = map(homodist_pool) do inst
            get_homogeneous_distribution_from_costs_sampled(
                inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
                max_num_edges=max_edges,
                rng=MersenneTwister(SEED + hash(inst.num_edges) + n_target))
        end
        sampled_avg = average_distributions(sampled_homodists)
        P_emp = compute_empirical_P(homodist_pool, max_edges)
        m_hd = size(sampled_avg, 1) - 1
        P_vals_samp = P_emp[1:min(m_hd+1, length(P_emp))]
        if length(P_vals_samp) < m_hd + 1
            P_vals_samp = vcat(P_vals_samp, zeros(m_hd + 1 - length(P_vals_samp)))
        end

        for p in P_VALUES
            println("    p=$p:")

            # Get transfer params
            transfer_γ, transfer_β = get_transfer_params(source_insts, p)

            # Get SampN+EmpP params
            samp_γ, samp_β, _ = optimize_via_proxy(sampled_avg, P_vals_samp, n_target, p)

            results_p = Dict{String, Vector{Float64}}()

            # === TRANSFER (no refine) ===
            transfer_ratios = map(enumerate(eval_insts)) do (i, inst)
                qaoa_expectation(inst.costs, n_target, transfer_γ, transfer_β) / eval_optimal[i]
            end
            results_p["Transfer"] = transfer_ratios
            @printf("      Transfer:       mean=%.4f  std=%.4f\n", mean(transfer_ratios), std(transfer_ratios))

            # === TRANSFER + REFINE ===
            t_tr_ref = @elapsed begin
                transfer_refine_ratios = map(enumerate(eval_insts)) do (i, inst)
                    rng_h = MersenneTwister(SEED + i * 3000 + hash(gt))
                    _, _, best_exp = hybrid_optimize(inst.costs, n_target, transfer_γ, transfer_β; rng=rng_h)
                    best_exp / eval_optimal[i]
                end
            end
            results_p["Transfer+Refine"] = transfer_refine_ratios
            @printf("      Transfer+Refine: mean=%.4f  std=%.4f  (%.1fs)\n",
                mean(transfer_refine_ratios), std(transfer_refine_ratios), t_tr_ref)

            # === SAMPN+EMPP (no refine) ===
            samp_ratios = map(enumerate(eval_insts)) do (i, inst)
                qaoa_expectation(inst.costs, n_target, samp_γ, samp_β) / eval_optimal[i]
            end
            results_p["SampN+EmpP"] = samp_ratios
            @printf("      SampN+EmpP:     mean=%.4f  std=%.4f\n", mean(samp_ratios), std(samp_ratios))

            # === SAMPN+EMPP + REFINE ===
            t_samp_ref = @elapsed begin
                samp_refine_ratios = map(enumerate(eval_insts)) do (i, inst)
                    rng_h = MersenneTwister(SEED + i * 4000 + hash(gt))
                    _, _, best_exp = hybrid_optimize(inst.costs, n_target, samp_γ, samp_β; rng=rng_h)
                    best_exp / eval_optimal[i]
                end
            end
            results_p["SampN+Refine"] = samp_refine_ratios
            @printf("      SampN+Refine:   mean=%.4f  std=%.4f  (%.1fs)\n",
                mean(samp_refine_ratios), std(samp_refine_ratios), t_samp_ref)

            # === RANDOM + REFINE ===
            t_rand = @elapsed begin
                random_refine_ratios = map(enumerate(eval_insts)) do (i, inst)
                    rng_r = MersenneTwister(SEED + i * 5000 + hash(gt))
                    best_exp = -Inf
                    for _ in 1:RANDOM_RESTARTS
                        init_γ = rand(rng_r, p) .* 1.6
                        init_β = rand(rng_r, p) .* (π/2)
                        _, _, ref_exp = refine_qaoa_from(inst.costs, n_target, init_γ, init_β)
                        best_exp = max(best_exp, ref_exp)
                    end
                    best_exp / eval_optimal[i]
                end
            end
            results_p["Random+Refine"] = random_refine_ratios
            @printf("      Random+Refine:  mean=%.4f  std=%.4f  (%.1fs)\n",
                mean(random_refine_ratios), std(random_refine_ratios), t_rand)

            all_results[(gt, n_target, p)] = results_p
        end
    end
end

#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY — Mean Approximation Ratios")
println("=" ^ 80)

for gt in graph_types
    println("\n--- $gt ---")
    @printf("%-4s %-4s %-10s %-10s %-10s %-10s %-10s %-12s %-12s\n",
        "n", "p", "Transfer", "Tr+Ref", "SampN", "Samp+Ref", "Rnd+Ref", "Gap(SR-TR)", "Gap(SR-T)")
    println("-" ^ 95)
    for n in N_TARGETS, p in P_VALUES
        key = (gt, n, p)
        haskey(all_results, key) || continue
        r = all_results[key]
        t = mean(r["Transfer"])
        tr = mean(r["Transfer+Refine"])
        s = mean(r["SampN+EmpP"])
        sr = mean(r["SampN+Refine"])
        rr = mean(r["Random+Refine"])
        gap_sr_tr = sr - tr
        gap_sr_t = sr - t
        marker = gap_sr_tr > 0.005 ? " ***" : gap_sr_tr > 0 ? " *" : ""
        @printf("%-4d %-4d %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f %+.4f     %+.4f%s\n",
            n, p, t, tr, s, sr, rr, gap_sr_tr, gap_sr_t, marker)
    end
end

println("\nDone!")
