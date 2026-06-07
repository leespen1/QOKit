#=
transfer_refine_investigation.jl — Compare warmstart sources for hybrid refinement.

Critical question: Is the proxy warmstart necessary, or does Transfer+Refine
achieve the same result? If Transfer+Refine ≈ Proxy+Refine, the proxy step
is unnecessary overhead.

Methods:
  1. Transfer (no refinement)
  2. PaperProxy (no refinement)
  3. Transfer+Refine (transfer params → coord descent on target)
  4. Proxy+Refine (proxy params → coord descent on target)
  5. Random+Refine (random start → coord descent)

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
const GRID_SIZE_P1 = 50
const GRID_SIZE_RAMP = 10
const REFINE_RESTARTS = 3
const REFINE_PERTURBATION = 0.2
const RANDOM_RESTARTS = 10
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

function refine_qaoa_from(costs, n, γs_init, βs_init)
    cγ, cβ = copy(γs_init), copy(βs_init); p = length(γs_init)
    ce = qaoa_expectation(costs, n, cγ, cβ)
    for ss in [0.3, 0.15, 0.07, 0.03, 0.01]
        for pi in 1:(2p), d in [-2,-1,-0.5,0.5,1,2] .* ss
            tγ, tβ = copy(cγ), copy(cβ)
            pi <= p ? (tγ[pi] = max(0, tγ[pi] + d)) : (tβ[pi-p] = clamp(tβ[pi-p] + d, 0, π/2))
            te = qaoa_expectation(costs, n, tγ, tβ)
            if te > ce; ce = te; cγ, cβ = tγ, tβ; end
        end
    end
    return cγ, cβ, ce
end

function hybrid_optimize(costs, n, init_γ, init_β; n_restarts=REFINE_RESTARTS,
    perturb=REFINE_PERTURBATION, rng=Random.default_rng())
    bγ, bβ, be = refine_qaoa_from(costs, n, init_γ, init_β)
    for _ in 1:n_restarts
        tγ = init_γ .+ randn(rng, length(init_γ)) .* perturb
        tβ = init_β .+ randn(rng, length(init_β)) .* perturb
        tγ = max.(tγ, 0); tβ = clamp.(tβ, 0, π/2)
        rγ, rβ, re = refine_qaoa_from(costs, n, tγ, tβ)
        if re > be; be = re; bγ, bβ = rγ, rβ; end
    end
    return bγ, bβ, be
end

function optimize_via_proxy(homodist, P_vals, n, p)
    if p == 1
        K = GRID_SIZE_P1^2; γm = zeros(K, 1); βm = zeros(K, 1); idx = 0
        for γ in range(0.02, 2.0, length=GRID_SIZE_P1),
            β in range(0.01, π/2-0.01, length=GRID_SIZE_P1)
            idx += 1; γm[idx, 1] = γ/π; βm[idx, 1] = β/π
        end
        Qs = QAOA_proxy_multi(homodist, γm, βm; pi_units=true)
        exps = vec(expectation(Qs[end], P_vals, n)); bi = argmax(exps)
        return [γm[bi,1]*π], [βm[bi,1]*π], exps[bi]
    else
        gs = GRID_SIZE_RAMP; K = gs^4; γm = zeros(K, p); βm = zeros(K, p); idx = 0
        for γ₁ in range(0.02,0.40,length=gs), γf in range(0.10,0.70,length=gs),
            β₁ in range(0.05,0.45,length=gs), βf in range(0.01,0.25,length=gs)
            idx += 1; γs, βs = linear_ramp(γ₁, γf, β₁, βf, p)
            γm[idx, :] .= γs; βm[idx, :] .= βs
        end
        Qs = QAOA_proxy_multi(homodist, γm, βm; pi_units=true)
        exps = vec(expectation(Qs[end], P_vals, n)); bi = argmax(exps)
        return γm[bi,:] .* π, βm[bi,:] .* π, exps[bi]
    end
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
println("TRANSFER+REFINE vs PROXY+REFINE")
println("=" ^ 80)

graph_types = ["ER", "BA", "WS"]
methods = ["Transfer", "PaperProxy", "Xfer+Refine", "Proxy+Refine", "Random+Refine"]

source = Dict(gt => generate_instances(gt, N_SMALL, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 1)) for gt in graph_types)
target = Dict(gt => generate_instances(gt, N_LARGE, NUM_INSTANCES;
    rng=MersenneTwister(SEED + hash(gt) + 2)) for gt in graph_types)
tgt_opt = Dict(gt => [maxcut_optimal(inst.costs) for inst in target[gt]] for gt in graph_types)

results = Dict(gt => Dict(m => Dict{Int, Vector{Float64}}() for m in methods) for gt in graph_types)

for gt in graph_types
    println("\n--- $gt ---")
    for p in P_VALUES
        println("  p=$p:")

        # Transfer params (same for Transfer and Xfer+Refine)
        if p == 1
            rng_opt = MersenneTwister(SEED + hash(gt) + p * 100)
            sp = map(source[gt]) do inst
                γs, βs, _ = optimize_qaoa_cd(inst.costs, N_SMALL, 1, N_RESTARTS; rng=rng_opt)
                (γs=γs, βs=βs)
            end
            med_γ = [median([s.γs[1] for s in sp])]
            med_β = [median([s.βs[1] for s in sp])]
        else
            sr = map(source[gt]) do inst
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

        # Transfer
        tr = map(enumerate(target[gt])) do (i, inst)
            qaoa_expectation(inst.costs, N_LARGE, med_γ, med_β) / tgt_opt[gt][i]
        end
        results[gt]["Transfer"][p] = tr

        # PaperProxy
        pp = map(enumerate(target[gt])) do (i, inst)
            pe = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, N_LARGE)
            proxy = PaperProxy(inst.num_edges, N_LARGE, pe)
            hd = cpu_compute_homodist(proxy)
            pv = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            bγ, bβ, _ = optimize_via_proxy(hd, pv, N_LARGE, p)
            qaoa_expectation(inst.costs, N_LARGE, bγ, bβ) / tgt_opt[gt][i]
        end
        results[gt]["PaperProxy"][p] = pp

        # Transfer + Refine
        xr = map(enumerate(target[gt])) do (i, inst)
            rng_h = MersenneTwister(SEED + i * 1000)
            _, _, be = hybrid_optimize(inst.costs, N_LARGE, med_γ, med_β; rng=rng_h)
            be / tgt_opt[gt][i]
        end
        results[gt]["Xfer+Refine"][p] = xr

        # Proxy + Refine
        pr = map(enumerate(target[gt])) do (i, inst)
            pe = gt == "ER" ? ER_P_EDGE : effective_edge_probability(inst.num_edges, N_LARGE)
            proxy = PaperProxy(inst.num_edges, N_LARGE, pe)
            hd = cpu_compute_homodist(proxy)
            pv = [P_cost_distribution(proxy, c) for c in 0:inst.num_edges]
            pγ, pβ, _ = optimize_via_proxy(hd, pv, N_LARGE, p)
            rng_h = MersenneTwister(SEED + i * 1000)
            _, _, be = hybrid_optimize(inst.costs, N_LARGE, pγ, pβ; rng=rng_h)
            be / tgt_opt[gt][i]
        end
        results[gt]["Proxy+Refine"][p] = pr

        # Random + Refine
        rr = map(enumerate(target[gt])) do (i, inst)
            rng_r = MersenneTwister(SEED + i * 2000); be = -Inf
            for _ in 1:RANDOM_RESTARTS
                iγ = rand(rng_r, p) .* 1.6; iβ = rand(rng_r, p) .* (π/2)
                _, _, re = refine_qaoa_from(inst.costs, N_LARGE, iγ, iβ)
                be = max(be, re)
            end
            be / tgt_opt[gt][i]
        end
        results[gt]["Random+Refine"][p] = rr

        for m in methods
            @printf("    %-15s mean=%.4f  std=%.4f\n", m, mean(results[gt][m][p]), std(results[gt][m][p]))
        end
    end
end

#==============================================================================#
#   PLOTS                                                                      #
#==============================================================================#

println("\n--- Plotting ---")
mc = Dict("Transfer"=>:steelblue, "PaperProxy"=>:coral, "Xfer+Refine"=>:orange,
    "Proxy+Refine"=>:mediumseagreen, "Random+Refine"=>:gray60)

fig = Figure(size=(500 * length(graph_types), 550))
for (gi, gt) in enumerate(graph_types)
    ax = Axis(fig[1, gi], xlabel="QAOA Depth p", ylabel="Approx Ratio",
        title="$gt (n=$N_SMALL→$N_LARGE)", xticks=P_VALUES)
    for p in P_VALUES
        nm = length(methods); tw = 0.85; sw = tw / nm
        for (mi, m) in enumerate(methods)
            offset = (mi - (nm+1)/2) * sw
            vals = results[gt][m][p]
            boxplot!(ax, fill(Float64(p)+offset, length(vals)), vals,
                color=mc[m], width=sw*0.8, label=(p==P_VALUES[1] ? m : nothing))
        end
    end
    if gi == 1; axislegend(ax, position=:rb, labelsize=8); end
end
Label(fig[0, :], "Warmstart Comparison: Transfer vs Proxy\n$NUM_INSTANCES instances",
    fontsize=14, font=:bold)
save_figure(fig, "transfer_vs_proxy_refine.png")

#==============================================================================#
#   SUMMARY                                                                    #
#==============================================================================#

println("\n" * "=" ^ 80)
println("KEY COMPARISON: Xfer+Refine vs Proxy+Refine")
println("=" ^ 80)
for p in P_VALUES
    @printf("  p=%d:", p)
    for gt in graph_types
        xr = mean(results[gt]["Xfer+Refine"][p])
        pr = mean(results[gt]["Proxy+Refine"][p])
        @printf("  %s: Xfer+R=%.4f Proxy+R=%.4f diff=%+.4f", gt, xr, pr, pr-xr)
    end
    println()
end

println("\nAll methods at p=5:")
@printf("  %-15s", "Method")
for gt in graph_types; @printf("  %-12s", gt); end
println()
for m in methods
    @printf("  %-15s", m)
    for gt in graph_types
        @printf("  %.4f±%.4f", mean(results[gt][m][5]), std(results[gt][m][5]))
    end
    println()
end

println("\nDone!")
