#=
instance_count_sparse_er.jl — How many graph instances does SampN+EmpP need
on sparse ER(p=0.2, 0.3)?

Sparse ER has fewer edges per instance, so each instance provides less
information. We need to verify that 10-20 instances is sufficient.

Test at n=18 with varying instance counts: 1, 3, 5, 10, 20, 40.

Started: 2026-04-02
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_TARGET = 18
const P_EDGES = [0.2, 0.3, 0.5]
const INSTANCE_COUNTS = [1, 3, 5, 10, 20, 40]
const NUM_EVAL = 5
const SAMPLES_PER_COST = 20
const P_VALUES = [1, 3]
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const SEED = 42
const NUM_REPEATS = 5  # repeat with different random seeds

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
    if p == 1
        K = GRID_SIZE_P1^2
        gm = zeros(K, 1); bm = zeros(K, 1); idx = 0
        for g in range(0.02, 2.0, length=GRID_SIZE_P1), b in range(0.01, pi/2-0.01, length=GRID_SIZE_P1)
            idx += 1; gm[idx,1] = g/pi; bm[idx,1] = b/pi
        end
        Qs = QAOA_proxy_multi(homodist, gm, bm; pi_units=true)
        exps = vec(expectation(Qs[end], P_vals, n))
        bi = argmax(exps)
        return [gm[bi,1]*pi], [bm[bi,1]*pi]
    else
        gs = GRID_SIZE_RAMP; K = gs^4
        gm = zeros(K, p); bm = zeros(K, p); idx = 0
        for g1 in range(0.02,0.40,length=gs), gf in range(0.10,0.70,length=gs),
            b1 in range(0.05,0.45,length=gs), bf in range(0.01,0.25,length=gs)
            idx += 1
            gs_v, bs_v = linear_ramp(g1, gf, b1, bf, p)
            gm[idx,:] .= gs_v; bm[idx,:] .= bs_v
        end
        Qs = QAOA_proxy_multi(homodist, gm, bm; pi_units=true)
        exps = vec(expectation(Qs[end], P_vals, n))
        bi = argmax(exps)
        return gm[bi,:] .* pi, bm[bi,:] .* pi
    end
end

#==============================================================================#
#   MAIN                                                                       #
#==============================================================================#

println("=" ^ 80)
println("Instance Count Sensitivity on Sparse ER (n=$N_TARGET)")
println("=" ^ 80)

# Generate eval instances (fixed)
eval_insts = Dict{Float64, Vector}()
eval_opts = Dict{Float64, Vector{Float64}}()
for pe in P_EDGES
    rng = MersenneTwister(SEED + hash(pe) + N_TARGET * 100)
    eval_insts[pe] = [generate_er_instance(N_TARGET, pe; rng) for _ in 1:NUM_EVAL]
    eval_opts[pe] = [maxcut_optimal(inst.costs) for inst in eval_insts[pe]]
end

for pe in P_EDGES
    println("\n--- p_edge=$pe ---")
    for p in P_VALUES
        println("  p=$p:")
        @printf("    %-12s", "Instances")
        for r in 1:NUM_REPEATS
            @printf("  seed%-2d  ", r)
        end
        @printf("  mean     std\n")

        for n_inst in INSTANCE_COUNTS
            @printf("    %-12d", n_inst)
            ratios_per_repeat = Float64[]

            for rep in 1:NUM_REPEATS
                rng = MersenneTwister(SEED + hash(pe) + N_TARGET * 200 + rep * 1000)
                hd_insts = [generate_er_instance(N_TARGET, pe; rng) for _ in 1:n_inst]
                max_edges = maximum(inst.num_edges for inst in hd_insts)

                sampled_hds = map(hd_insts) do inst
                    get_homogeneous_distribution_from_costs_sampled(
                        inst.costs, inst.num_edges, inst.num_vertices, SAMPLES_PER_COST;
                        max_num_edges=max_edges,
                        rng=MersenneTwister(SEED + hash(pe) + hash(inst.num_edges) + N_TARGET + rep))
                end
                avg_hd = average_distributions(sampled_hds)
                P_emp = compute_empirical_P(hd_insts, max_edges)
                m_hd = size(avg_hd, 1) - 1
                Pv = P_emp[1:min(m_hd+1, length(P_emp))]
                if length(Pv) < m_hd + 1; Pv = vcat(Pv, zeros(m_hd+1-length(Pv))); end

                bg, bb = optimize_via_proxy(avg_hd, Pv, N_TARGET, p)
                mean_ratio = mean(
                    qaoa_expectation(inst.costs, N_TARGET, bg, bb) / eval_opts[pe][i]
                    for (i, inst) in enumerate(eval_insts[pe]))
                push!(ratios_per_repeat, mean_ratio)
                @printf("  %.4f  ", mean_ratio)
            end
            @printf("  %.4f  %.4f\n", mean(ratios_per_repeat), std(ratios_per_repeat))
        end
    end
end

println("\nDone!")
