#=
paperproxy_advantage_investigation.jl — Why does PaperProxy+p_eff outperform
empirical homodist on non-ER graphs at high depth?

Hypotheses:
  H1: PaperProxy landscape is smoother, leading to better optimization
  H2: p_eff accidentally matches BA/WS cost structure well
  H3: The binomial P(c') is a better regularizer than empirical P
  H4: Analytical formula provides implicit instance-averaging that the
      empirical approach (even multi-instance) doesn't capture fully

Tests:
  1. Compare proxy landscapes (heatmaps at p=1)
  2. Vary p_eff for PaperProxy on non-ER — is there an optimal p_eff?
  3. PaperProxy N + Empirical P vs PaperProxy N + Binomial P
  4. Increase number of homodist instances (5 → 100) for EmpN+EmpP

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf
using Distributions: Binomial, pdf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_TARGET = 12
const NUM_EVAL_INSTANCES = 20
const N_SOURCE = 9
const SEED = 42

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

const GRID_SIZE_P1 = 60
const GRID_SIZE_RAMP = 10

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

function optimize_via_proxy_p1(homodist, P_vals, n)
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
    return [γ_matrix[best_idx, 1] * π], [β_matrix[best_idx, 1] * π], exps[best_idx], reshape(exps, GRID_SIZE_P1, GRID_SIZE_P1), γ_range, β_range
end

function optimize_via_proxy_ramp(homodist, P_vals, n, p)
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

function print_stats(label, values)
    @printf("  %-55s  mean=%.4f  std=%.4f\n", label, mean(values), std(values))
end


#==============================================================================#
#   TEST 1: Vary p_eff for PaperProxy on non-ER                               #
#==============================================================================#

println("=" ^ 80)
println("TEST 1: PaperProxy with varying p_eff on non-ER graphs")
println("=" ^ 80)

graph_types = ["BA", "WS"]

eval_instances = Dict(gt => generate_instances(gt, N_TARGET, NUM_EVAL_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 2)) for gt in graph_types)
eval_optimal = Dict(gt => [maxcut_optimal(inst.costs) for inst in eval_instances[gt]]
    for gt in graph_types)

p_eff_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

for gt in graph_types
    println("\n--- $gt ---")
    for p in [1, 3]
        println("  p=$p:")
        for p_eff in p_eff_values
            ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                proxy = PaperProxy(inst.num_edges, N_TARGET, p_eff)
                homodist = cpu_compute_homodist(proxy)
                P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
                if p == 1
                    best_γ, best_β, _, _, _, _ = optimize_via_proxy_p1(homodist, P_vals, N_TARGET)
                else
                    best_γ, best_β, _ = optimize_via_proxy_ramp(homodist, P_vals, N_TARGET, p)
                end
                qaoa_expectation(inst.costs, N_TARGET, best_γ, best_β) / eval_optimal[gt][i]
            end
            @printf("    p_eff=%.1f: mean=%.4f  std=%.4f\n", p_eff, mean(ratios), std(ratios))
        end
        # Also show natural p_eff
        nat_peff = mean(effective_edge_probability(inst.num_edges, N_TARGET) for inst in eval_instances[gt])
        @printf("    (natural p_eff ≈ %.3f)\n", nat_peff)
    end
end


#==============================================================================#
#   TEST 2: EmpN+EmpP with varying number of instances                        #
#==============================================================================#

println("\n" * "=" ^ 80)
println("TEST 2: EmpN+EmpP with varying number of homodist instances")
println("=" ^ 80)

instance_counts = [5, 10, 20, 50, 100]

for gt in graph_types
    println("\n--- $gt ---")

    # Generate a large pool of instances
    pool = generate_instances(gt, N_TARGET, 100;
        rng=MersenneTwister(SEED + hash(gt) + 500))
    max_edges_pool = maximum(inst.num_edges for inst in pool)

    for p in [1, 3]
        println("  p=$p:")
        for num_inst in instance_counts
            subset = pool[1:num_inst]
            homodists = map(subset) do inst
                get_homogeneous_distribution_from_costs_direct(
                    inst.costs, inst.num_edges, inst.num_vertices;
                    max_num_edges=max_edges_pool)
            end
            avg_hd = average_distributions(homodists)
            P_vals = compute_empirical_P(subset, max_edges_pool)
            m_hd = size(avg_hd, 1) - 1
            P_v = P_vals[1:min(m_hd+1, length(P_vals))]
            if length(P_v) < m_hd + 1
                P_v = vcat(P_v, zeros(m_hd + 1 - length(P_v)))
            end

            if p == 1
                best_γ, best_β, _, _, _, _ = optimize_via_proxy_p1(avg_hd, P_v, N_TARGET)
            else
                best_γ, best_β, _ = optimize_via_proxy_ramp(avg_hd, P_v, N_TARGET, p)
            end

            ratios = map(enumerate(eval_instances[gt])) do (i, inst)
                qaoa_expectation(inst.costs, N_TARGET, best_γ, best_β) / eval_optimal[gt][i]
            end
            @printf("    N_inst=%3d: mean=%.4f  std=%.4f\n", num_inst, mean(ratios), std(ratios))
        end
    end
end


#==============================================================================#
#   TEST 3: PaperProxy N + Empirical P vs PaperProxy N + Binomial P           #
#==============================================================================#

println("\n" * "=" ^ 80)
println("TEST 3: PaperProxy N with different P(c') sources")
println("=" ^ 80)

# Get empirical P from pool instances
emp_P_pool = Dict{String, Vector{Float64}}()
for gt in graph_types
    pool = generate_instances(gt, N_TARGET, 50;
        rng=MersenneTwister(SEED + hash(gt) + 600))
    max_edges_pool = maximum(inst.num_edges for inst in pool)
    emp_P_pool[gt] = compute_empirical_P(pool, max_edges_pool)
end

for gt in graph_types
    println("\n--- $gt ---")
    for p in [1, 3]
        println("  p=$p:")

        # PaperProxy N + Binomial P (standard)
        ratios_binomial = map(enumerate(eval_instances[gt])) do (i, inst)
            p_eff = effective_edge_probability(inst.num_edges, N_TARGET)
            proxy = PaperProxy(inst.num_edges, N_TARGET, p_eff)
            homodist = cpu_compute_homodist(proxy)
            P_vals = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            if p == 1
                best_γ, best_β, _, _, _, _ = optimize_via_proxy_p1(homodist, P_vals, N_TARGET)
            else
                best_γ, best_β, _ = optimize_via_proxy_ramp(homodist, P_vals, N_TARGET, p)
            end
            qaoa_expectation(inst.costs, N_TARGET, best_γ, best_β) / eval_optimal[gt][i]
        end
        print_stats("PaperN + BinomialP", ratios_binomial)

        # PaperProxy N + Empirical P
        ratios_emp = map(enumerate(eval_instances[gt])) do (i, inst)
            p_eff = effective_edge_probability(inst.num_edges, N_TARGET)
            proxy = PaperProxy(inst.num_edges, N_TARGET, p_eff)
            homodist = cpu_compute_homodist(proxy)
            P_full = emp_P_pool[gt]
            m_proxy = inst.num_edges
            P_vals = P_full[1:min(m_proxy+1, length(P_full))]
            if length(P_vals) < m_proxy + 1
                P_vals = vcat(P_vals, zeros(m_proxy + 1 - length(P_vals)))
            end
            if p == 1
                best_γ, best_β, _, _, _, _ = optimize_via_proxy_p1(homodist, P_vals, N_TARGET)
            else
                best_γ, best_β, _ = optimize_via_proxy_ramp(homodist, P_vals, N_TARGET, p)
            end
            qaoa_expectation(inst.costs, N_TARGET, best_γ, best_β) / eval_optimal[gt][i]
        end
        print_stats("PaperN + EmpiricalP", ratios_emp)
    end
end


#==============================================================================#
#   TEST 4: Landscape Comparison (p=1 heatmaps)                               #
#==============================================================================#

println("\n" * "=" ^ 80)
println("TEST 4: Proxy Landscape Comparison (p=1)")
println("=" ^ 80)

# Pick one representative instance per graph type
fig = Figure(size=(400 * 3, 400 * length(graph_types)))

for (gi, gt) in enumerate(graph_types)
    inst = eval_instances[gt][1]
    p_eff = effective_edge_probability(inst.num_edges, N_TARGET)

    # True QAOA landscape
    γ_range = range(0.02, 2.0, length=GRID_SIZE_P1)
    β_range = range(0.01, π/2 - 0.01, length=GRID_SIZE_P1)
    true_landscape = [qaoa_expectation(inst.costs, N_TARGET, [γ], [β]) / maxcut_optimal(inst.costs)
        for γ in γ_range, β in β_range]

    # PaperProxy landscape
    proxy = PaperProxy(inst.num_edges, N_TARGET, p_eff)
    homodist_pp = cpu_compute_homodist(proxy)
    P_pp = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
    _, _, _, pp_landscape, _, _ = optimize_via_proxy_p1(homodist_pp, P_pp, N_TARGET)

    # EmpN+EmpP landscape
    pool = generate_instances(gt, N_TARGET, 50; rng=MersenneTwister(SEED + hash(gt) + 700))
    max_edges_pool = maximum(i.num_edges for i in pool)
    homodists_emp = map(pool) do i
        get_homogeneous_distribution_from_costs_direct(
            i.costs, i.num_edges, i.num_vertices; max_num_edges=max_edges_pool)
    end
    avg_hd = average_distributions(homodists_emp)
    P_emp = compute_empirical_P(pool, max_edges_pool)
    m_hd = size(avg_hd, 1) - 1
    P_v = P_emp[1:min(m_hd+1, length(P_emp))]
    if length(P_v) < m_hd + 1
        P_v = vcat(P_v, zeros(m_hd + 1 - length(P_v)))
    end
    _, _, _, emp_landscape, _, _ = optimize_via_proxy_p1(avg_hd, P_v, N_TARGET)

    # Normalize proxy landscapes to same scale as true
    pp_normalized = pp_landscape ./ maximum(pp_landscape) .* maximum(true_landscape)
    emp_normalized = emp_landscape ./ maximum(emp_landscape) .* maximum(true_landscape)

    ax1 = Axis(fig[gi, 1], xlabel="γ", ylabel="β", title="$gt: True QAOA")
    heatmap!(ax1, collect(γ_range), collect(β_range), true_landscape)

    ax2 = Axis(fig[gi, 2], xlabel="γ", ylabel="β", title="$gt: PaperProxy (p_eff)")
    heatmap!(ax2, collect(γ_range), collect(β_range), pp_normalized)

    ax3 = Axis(fig[gi, 3], xlabel="γ", ylabel="β", title="$gt: EmpN+EmpP")
    heatmap!(ax3, collect(γ_range), collect(β_range), emp_normalized)
end

Label(fig[0, :], "Proxy Landscape Comparison (p=1, n=$N_TARGET)", fontsize=14, font=:bold)
save_figure(fig, "paperproxy_advantage_landscapes.png")


#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

println("""
Key findings:
1. TEST 1: Does PaperProxy performance vary with p_eff? If so, does the natural
   p_eff happen to be near-optimal?
2. TEST 2: Does increasing homodist instances help EmpN+EmpP catch up to PaperProxy?
3. TEST 3: Does mixing PaperProxy N with empirical P help or hurt? (Tests
   whether PaperProxy's advantage is in N or P)
4. TEST 4: Visual comparison of landscape shapes — is PaperProxy's landscape
   closer to the true landscape?
""")

println("Done!")
