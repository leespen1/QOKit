#=
E032 — AUDIT: is the heavy-tail "proxy wins" niche an artifact of a
handicapped (mean-rescaled) transfer baseline?

E022 rescaled transfer's normalized γ̂-grid by each instance's MEAN edge
weight. Under Pareto(1.5) the mean is tail-dominated, so that transfer is
crippled exactly where the proxy is claimed to win. This experiment repeats
the headline Pareto ER(0.5) n=16 cell (same seeds as E022) and compares the
binned proxy against transfer rescaled three ways:
  transfer_mean    — γ = γ̂ / mean(w)          (E022's baseline)
  transfer_median  — γ = γ̂ / median(w)         (robust location)
  transfer_coststd — γ = γ̂ / std(costs)        (robust cost-scale;
                     what the scale-free proxy implicitly uses)
The source picks the best NORMALIZED index γ̂ on its own landscape under the
SAME rescaling, so each is a fair mean/median/std-rescaled angle transfer.
If median/coststd transfer ≈ the proxy, the niche is a baseline artifact.

p=1 only, 10 instances (E022's seeds). Local: JULIA_NUM_THREADS=8 julia --project run.jl
=#

using JuliaQAOA
using Random: MersenneTwister, rand
using Base.Threads: @threads
using Statistics: mean, median, std

const SEED = 20260708          # E022's target seeds
const SRC_SEED = 20260709      # E022's source seeds
const SOURCE_N = 10
const N_SOURCES = 3
const K = 64
const S_N = 10
const NINST = 10
const N = 16
const GRID_LEN = 40

const P1_γ̂ = collect(range(0.0, 2π; length=GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=GRID_LEN))
pareto15(rng, m) = (1 .- rand(rng, m)).^(-1 / 1.5)

function weighted_costs(n, edges, w)
    costs = zeros(Float64, 1 << n)
    @threads for x in 0:((1 << n) - 1)
        c = 0.0
        for (k, (i, j)) in enumerate(edges)
            ((x >> i) & 1) != ((x >> j) & 1) && (c += w[k])
        end
        costs[x + 1] = c
    end
    costs
end

function quantile_bins(costs, K)
    Nn = length(costs); order = sortperm(costs)
    labels = Vector{Int}(undef, Nn)
    for (rank, idx) in enumerate(order)
        labels[idx] = min(K - 1, div((rank - 1) * K, Nn))
    end
    sums = zeros(K); sizes = zeros(Int, K)
    for i in 1:Nn
        b = labels[i] + 1; sums[b] += costs[i]; sizes[b] += 1
    end
    labels, sums ./ sizes, sizes
end

function binned_proxy_argmax(Nd, binmeans, binsizes, n, scheds)
    Kk = length(binmeans)
    βset = Dict{Float64, Vector{ComplexF64}}()
    for (γs, βs) in scheds, β in βs
        haskey(βset, β) || (βset[β] = vec(get_β_factors([β], n)))
    end
    vals = zeros(length(scheds))
    @threads for k in eachindex(scheds)
        γs, βs = scheds[k]
        Q = fill(ComplexF64(1 / sqrt(2.0^n)), Kk)
        Qsrc = zeros(ComplexF64, Kk); Qnew = zeros(ComplexF64, Kk)
        for ℓ in eachindex(γs)
            f = βset[βs[ℓ]]; γ = γs[ℓ]
            @inbounds for b in 1:Kk
                Qsrc[b] = cis(-γ * binmeans[b] / 2) * Q[b]
            end
            fill!(Qnew, zero(ComplexF64))
            @inbounds for b in 1:Kk, d in 0:n
                w0 = f[d + 1] * Qsrc[b]
                for b′ in 1:Kk
                    Qnew[b′] += w0 * Nd[b′, d + 1, b]
                end
            end
            Q, Qnew = Qnew, Q
        end
        vals[k] = sum(binsizes[b] * abs2(Q[b]) * binmeans[b] for b in 1:Kk)
    end
    argmax(vals)
end

true_vals(costs, n, scheds) = begin
    v = zeros(length(scheds))
    @threads for k in eachindex(scheds)
        v[k] = qaoa_expectation(costs, n, scheds[k]...)
    end
    v
end

"Schedules with γ = γ̂ / scale over the shared normalized grid."
scheds_for(scale) = vec([([γ̂ / scale], [β]) for γ̂ in P1_γ̂, β in P1_β])

# scale estimators, keyed by name; take (weights, costs)
const SCALES = ["mean", "median", "coststd"]
scale_of(name, w, costs) = name == "mean" ? mean(w) :
                           name == "median" ? median(w) : std(costs)

function main()
    # sources: one per rep, store per-scale best normalized index
    src_idx = Dict(s => Int[] for s in SCALES)
    for rep in 1:N_SOURCES
        srng = MersenneTwister(SRC_SEED + 1000 + rep)   # ER(0.5) => family index 1
        sedges = erdos_renyi_edges(SOURCE_N, 0.5; rng=srng)
        sw = pareto15(srng, length(sedges))
        scosts = weighted_costs(SOURCE_N, sedges, sw)
        copt = maximum(scosts)
        for s in SCALES
            sv = true_vals(scosts, SOURCE_N, scheds_for(scale_of(s, sw, scosts))) ./ copt
            push!(src_idx[s], argmax(sv))
        end
    end

    rows = String[]
    for inst in 1:NINST
        seed = SEED + 100 * N + inst
        rng = MersenneTwister(seed)
        edges = erdos_renyi_edges(N, 0.5; rng=rng)
        m = length(edges)
        w = pareto15(rng, m)
        costs = weighted_costs(N, edges, w)
        c_opt = maximum(costs)

        # proxy on a scale-free reference grid (mean, as E022; scale cancels in argmax over bins)
        pscheds = scheds_for(mean(w))
        tv_p = true_vals(costs, N, pscheds)
        ceil_ar = maximum(tv_p) / c_opt
        labels, binmeans, binsizes = quantile_bins(costs, K)
        N_samp = sampled_homogeneous_distribution(Float64.(labels), K - 1, N;
                     samples_per_class=S_N, rng=MersenneTwister(seed + 777))
        pidx = binned_proxy_argmax(N_samp, binmeans, binsizes, N, pscheds)
        push!(rows, join(Any[inst, m, "binned_sampled", tv_p[pidx] / c_opt,
                             ceil_ar - tv_p[pidx] / c_opt], ","))

        for s in SCALES
            sc = scale_of(s, w, costs)
            tsched = scheds_for(sc)
            tv = true_vals(costs, N, tsched)
            cap = maximum(tv) / c_opt   # this scaling's own ceiling (same 40x40 grid)
            for r in 1:N_SOURCES
                ar = tv[src_idx[s][r]] / c_opt
                push!(rows, join(Any[inst, m, "transfer_$(s)_$r", ar, cap - ar], ","))
            end
        end
        println("inst=$inst m=$m ceil=$(round(ceil_ar; digits=4))"); flush(stdout)
    end

    open(joinpath(@__DIR__, "results.csv"), "w") do io
        println(io, "instance,m,method,ar,regret")
        foreach(r -> println(io, r), rows)
    end
    println("E032 complete")
end

main()
