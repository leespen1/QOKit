#=
cross_type_transfer_depth.jl — Cross-type parameter transfer at p=1,2,3,5.

The smallworld investigation showed cross-type transfer works at p=1 (within
0.005 penalty). This tests whether it holds at higher depths with linear ramp.

Practical importance: If ER→BA transfer works, practitioners can optimize on
cheap ER source graphs and apply to any graph type.

Started: 2026-04-01
=#

include("smallworld_common.jl")
using Printf

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_SMALL = 9
const N_LARGE = 12
const NUM_INSTANCES = 20

const BA_M_ATTACH = 2
const WS_K = 4
const WS_P_REWIRE = 0.3
const ER_P_EDGE = 0.5

const P_VALUES = [1, 2, 3, 5]
const N_RESTARTS = 15
const GRID_SIZE_RAMP = 10
const SEED = 42

#==============================================================================#
#                    HELPERS                                                    #
#==============================================================================#

function optimize_qaoa_cd(costs, n, p, n_restarts; rng=Random.default_rng())
    best_exp = -Inf; best_γs = zeros(p); best_βs = zeros(p)
    for _ in 1:n_restarts
        γs = rand(rng, p) .* 1.6; βs = rand(rng, p) .* (π/2)
        cγ, cβ = copy(γs), copy(βs)
        ce = qaoa_expectation(costs, n, cγ, cβ)
        for ss in [0.3, 0.15, 0.07, 0.03, 0.01]
            for pi in 1:(2p), d in [-2,-1,-0.5,0.5,1,2] .* ss
                tγ, tβ = copy(cγ), copy(cβ)
                pi <= p ? (tγ[pi] = max(0, tγ[pi] + d)) : (tβ[pi-p] = clamp(tβ[pi-p] + d, 0, π/2))
                te = qaoa_expectation(costs, n, tγ, tβ)
                if te > ce; ce = te; cγ, cβ = tγ, tβ; end
            end
        end
        if ce > best_exp; best_exp = ce; best_γs, best_βs = copy(cγ), copy(cβ); end
    end
    return best_γs, best_βs, best_exp
end

function generate_instances(gt::String, n::Int, num::Int; rng=Random.default_rng())
    map(1:num) do _
        gt == "ER" ? generate_er_instance(n, ER_P_EDGE; rng) :
        gt == "BA" ? generate_ba_instance(n, BA_M_ATTACH; rng) :
        gt == "WS" ? generate_ws_instance(n, WS_K, WS_P_REWIRE; rng) :
        error("Unknown: $gt")
    end
end

#==============================================================================#
#   MAIN                                                                       #
#==============================================================================#

println("=" ^ 80)
println("CROSS-TYPE TRANSFER AT HIGHER DEPTHS")
println("Source n=$N_SMALL, Target n=$N_LARGE, $NUM_INSTANCES instances")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]

# Generate all instances
source = Dict(gt => generate_instances(gt, N_SMALL, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 1)) for gt in graph_types)
target = Dict(gt => generate_instances(gt, N_LARGE, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 2)) for gt in graph_types)
tgt_opt = Dict(gt => [maxcut_optimal(inst.costs) for inst in target[gt]] for gt in graph_types)

# results[src_type][tgt_type][p] = Vector{Float64}
results = Dict(s => Dict(t => Dict{Int, Vector{Float64}}()
    for t in graph_types) for s in graph_types)

for src in graph_types
    println("\n--- Source: $src ---")

    for p in P_VALUES
        # Compute source parameters
        if p == 1
            rng_opt = MersenneTwister(SEED + hash(src) + 100)
            sp = map(source[src]) do inst
                γs, βs, _ = optimize_qaoa_cd(inst.costs, N_SMALL, 1, N_RESTARTS; rng=rng_opt)
                (γs=γs, βs=βs)
            end
            med_γ = [median([s.γs[1] for s in sp])]
            med_β = [median([s.βs[1] for s in sp])]
        else
            sr = map(source[src]) do inst
                bp = (0.0,0.0,0.0,0.0); be = -Inf; gs = GRID_SIZE_RAMP
                for γ₁ in range(0.02,0.40,length=gs), γf in range(0.10,0.70,length=gs),
                    β₁ in range(0.05,0.45,length=gs), βf in range(0.01,0.25,length=gs)
                    γp, βp = linear_ramp(γ₁, γf, β₁, βf, p)
                    ev = qaoa_expectation(inst.costs, N_SMALL, γp.*π, βp.*π)
                    if ev > be; be = ev; bp = (γ₁,γf,β₁,βf); end
                end
                bp
            end
            mr = Tuple(median([s[i] for s in sr]) for i in 1:4)
            γt, βt = linear_ramp(mr..., p)
            med_γ = γt .* π; med_β = βt .* π
        end

        println("  p=$p:")
        for tgt in graph_types
            ratios = map(enumerate(target[tgt])) do (i, inst)
                qaoa_expectation(inst.costs, N_LARGE, med_γ, med_β) / tgt_opt[tgt][i]
            end
            results[src][tgt][p] = ratios
            @printf("    %s→%s: mean=%.4f std=%.4f\n", src, tgt, mean(ratios), std(ratios))
        end
    end
end

#==============================================================================#
#   PLOTS                                                                      #
#==============================================================================#

println("\n--- Plotting ---")

# Heatmap for each depth
fig = Figure(size=(350 * length(P_VALUES), 350))
for (pi, p) in enumerate(P_VALUES)
    mat = [mean(results[s][t][p]) for s in graph_types, t in graph_types]
    ax = Axis(fig[1, pi], xlabel="Target", ylabel="Source",
        title="p=$p", xticks=(1:3, graph_types), yticks=(1:3, graph_types))
    hm = heatmap!(ax, 1:3, 1:3, mat, colormap=:YlOrRd,
        colorrange=(minimum(mat) - 0.02, maximum(mat) + 0.02))
    for i in 1:3, j in 1:3
        text!(ax, j, i, text=@sprintf("%.3f", mat[i,j]),
            align=(:center, :center), fontsize=12, color=:black)
    end
    Colorbar(fig[1, pi+length(P_VALUES)], hm, label="ApxRatio", width=12)
end
Label(fig[0, :], "Cross-Type Transfer (n=$N_SMALL→$N_LARGE, $NUM_INSTANCES instances)",
    fontsize=14, font=:bold)
save_figure(fig, "cross_type_transfer_depth.png")

#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

for p in P_VALUES
    println("\n  p=$p:")
    @printf("  %-8s", "Src\\Tgt")
    for t in graph_types; @printf("  %-10s", t); end
    println()
    println("  " * "-" ^ 40)
    for s in graph_types
        @printf("  %-8s", s)
        for t in graph_types
            @printf("  %.4f    ", mean(results[s][t][p]))
        end
        println()
    end
end

# Cross-type penalty (same-type - cross-type)
println("\n\nCross-Type Penalty (same_type - cross_type, positive = penalty):")
for p in P_VALUES
    println("  p=$p:")
    for s in graph_types
        for t in graph_types
            if s == t; continue; end
            same = mean(results[t][t][p])
            cross = mean(results[s][t][p])
            @printf("    %s→%s: same=%.4f cross=%.4f penalty=%.4f\n",
                s, t, same, cross, same - cross)
        end
    end
end

println("\nDone!")
