#=
E031 — The composed recipe: cheap init + surrogate polish + one check.

The map says: init from transfer/universal where family concentration
holds, from the binned proxy under high dispersion; E030 says: polish on
the surrogate within the ramp family. This experiment tests the composed,
fully deployable recipe (no target evaluations except one confirmation)
in three regimes:
  task 1: unweighted ER(0.5), p=3    (transfer's home turf)
  task 2: Pareto(1.5) ER(0.5), p=3   (the proxy's niche)
  task 3: lognormal(2.5) 3-regular, p=20 (the hard corner)
Inits: transfer_1 (one source, no target evals) and universal (pooled).
Each polished on the binned proxy's ramp landscape; one true evaluation
per stage. Reported against the 8^4 grid ceiling.

Output: results_task<ID>.csv
(family,n,instance,seed,m,law,init,stage,ar,regret_vs_gridceil)
Submit:  cd research/experiments/031_composed-recipe && sbatch run.sb
Smoke:   E31_SMOKE=1 julia --project research/experiments/031_composed-recipe/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister, rand, randn
using Base.Threads: @threads
using Statistics: mean

const SMOKE = get(ENV, "E31_SMOKE", "0") == "1"

const SEED = 20260720            # E028's convention: same instances
const SRC_SEED = 20260721
const SOURCE_N = 10
const N_SOURCES = SMOKE ? 1 : 3
const K = 64
const S_N = 10
const RAMP_LEN = SMOKE ? 4 : 8

const RAMP_γ̂ = collect(range(0.1, 3.2; length=RAMP_LEN))
const RAMP_β = collect(range(0.05, 0.8; length=RAMP_LEN))

const TASKS = [("ER(0.5)", "unit", 3, 10), ("ER(0.5)", "pareto15", 3, 10),
               ("3-regular", "logn2.5", 20, 5)]
const N_TARGET = SMOKE ? 12 : 16
const GENS = Dict(
    "ER(0.5)"   => (rng, n) -> erdos_renyi_edges(n, 0.5; rng),
    "3-regular" => (rng, n) -> random_regular_edges(n, 3; rng),
)
const LAWS = Dict("pareto15" => (rng, m) -> (1 .- rand(rng, m)).^(-1/1.5),
                  "logn2.5"  => (rng, m) -> exp.(2.5 .* randn(rng, m)),
                  "unit"     => (rng, m) -> ones(m))

function weighted_costs(n, edges, w)
    costs = zeros(Float64, 1 << n)
    @threads for x in 0:((1 << n) - 1)
        c = 0.0
        for (k, (i, j)) in enumerate(edges)
            if ((x >> i) & 1) != ((x >> j) & 1)
                c += w[k]
            end
        end
        costs[x + 1] = c
    end
    costs
end

function quantile_bins(costs, K)
    N = length(costs)
    order = sortperm(costs)
    labels = Vector{Int}(undef, N)
    for (rank, idx) in enumerate(order)
        labels[idx] = min(K - 1, div((rank - 1) * K, N))
    end
    sums = zeros(K); sizes = zeros(Int, K)
    for i in 1:N
        b = labels[i] + 1
        sums[b] += costs[i]; sizes[b] += 1
    end
    @assert all(>(0), sizes)
    labels, sums ./ sizes, sizes
end

function binned_proxy_argmax(N, binmeans, binsizes, n, schedules)
    K = length(binmeans)
    βset = Dict{Float64, Vector{ComplexF64}}()
    for (γs, βs) in schedules, β in βs
        haskey(βset, β) || (βset[β] = vec(get_β_factors([β], n)))
    end
    vals = zeros(length(schedules))
    nt = Threads.nthreads()
    @threads :static for t in 1:nt
        Q = zeros(ComplexF64, K); Qsrc = zeros(ComplexF64, K); Qnew = zeros(ComplexF64, K)
        for k in t:nt:length(schedules)
            γs, βs = schedules[k]
            fill!(Q, ComplexF64(1 / sqrt(2.0^n)))
            for ℓ in eachindex(γs)
                f = βset[βs[ℓ]]; γ = γs[ℓ]
                @inbounds for b in 1:K
                    Qsrc[b] = cis(-γ * binmeans[b] / 2) * Q[b]
                end
                fill!(Qnew, zero(ComplexF64))
                @inbounds for b in 1:K
                    w0 = Qsrc[b]
                    for d in 0:n
                        w = f[d + 1] * w0
                        for b′ in 1:K
                            Qnew[b′] += w * N[b′, d + 1, b]
                        end
                    end
                end
                Q, Qnew = Qnew, Q
            end
            s = 0.0
            @inbounds for b in 1:K
                s += binsizes[b] * abs2(Q[b]) * binmeans[b]
            end
            vals[k] = s
        end
    end
    argmax(vals)
end

"Binned proxy predicted value for ONE schedule (allocation-light)."
function binned_proxy_value(N, binmeans, binsizes, n, γs, βs)
    K = length(binmeans)
    Q = fill(ComplexF64(1 / sqrt(2.0^n)), K)
    Qsrc = zeros(ComplexF64, K); Qnew = zeros(ComplexF64, K)
    for ℓ in eachindex(γs)
        f = vec(get_β_factors([βs[ℓ]], n)); γ = γs[ℓ]
        @inbounds for b in 1:K
            Qsrc[b] = cis(-γ * binmeans[b] / 2) * Q[b]
        end
        fill!(Qnew, zero(ComplexF64))
        @inbounds for b in 1:K
            w0 = Qsrc[b]
            for d in 0:n
                w = f[d + 1] * w0
                for b′ in 1:K
                    Qnew[b′] += w * N[b′, d + 1, b]
                end
            end
        end
        Q, Qnew = Qnew, Q
    end
    sum(binsizes[b] * abs2(Q[b]) * binmeans[b] for b in 1:K)
end

function true_vals(costs, n, scheds)
    vals = zeros(length(scheds))
    @threads for k in eachindex(scheds)
        γs, βs = scheds[k]
        vals[k] = qaoa_expectation(costs, n, γs, βs)
    end
    vals
end

"Compass (pattern) search maximizing f over θ, coordinate steps halving on failure."
function compass_search(f, θ0, step0; min_step=1e-3, max_evals=SMOKE ? 300 : 2000)
    θ = copy(θ0); fθ = f(θ)
    steps = copy(step0)
    evals = 1
    while evals < max_evals && maximum(steps) > min_step
        improved = false
        for i in eachindex(θ), s in (steps[i], -steps[i])
            evals >= max_evals && break
            θi = copy(θ); θi[i] += s
            fi = f(θi); evals += 1
            if fi > fθ
                θ, fθ = θi, fi
                improved = true
            end
        end
        improved || (steps .*= 0.5)
    end
    θ, fθ
end

function main()
    tid = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
    fam, law, P, NINST0 = SMOKE ? ("ER(0.5)", "pareto15", 5, 1) : TASKS[tid == 0 ? 1 : tid]
    NINST = NINST0
    n = N_TARGET
    drawlaw = LAWS[law]
    println("task $tid: $fam law=$law n=$n p=$P")

    scheds_for(w̄) = begin
        combos = vec([(g1 / w̄, gf / w̄, b1, bf) for g1 in RAMP_γ̂, gf in RAMP_γ̂,
                                                    b1 in RAMP_β, bf in RAMP_β])
        combos, [linear_ramp(c..., P) for c in combos]
    end

    # sources (E028 convention)
    srcv = Vector{Vector{Float64}}()
    src_combos = nothing
    for rep in 1:N_SOURCES
        srng = MersenneTwister(SRC_SEED + 1000 * (fam == "ER(0.5)" ? 1 : 2) + rep)
        sedges = GENS[fam](srng, SOURCE_N)
        sw = drawlaw(srng, length(sedges))
        scosts = weighted_costs(SOURCE_N, sedges, sw)
        combos, sscheds = scheds_for(mean(sw))
        src_combos = combos
        push!(srcv, true_vals(scosts, SOURCE_N, sscheds) ./ maximum(scosts))
    end
    transfer_idx = [argmax(v) for v in srcv]
    universal_idx = argmax(reduce(+, srcv))

    rows = String[]
    for inst in 1:NINST
        seed = SEED + 100 * n + inst
        rng = MersenneTwister(seed)
        edges = GENS[fam](rng, n)
        m = length(edges)
        w = drawlaw(rng, m)
        costs = weighted_costs(n, edges, w)
        c_opt = maximum(costs)
        w̄ = mean(w)
        combos, scheds = scheds_for(w̄)

        labels, binmeans, binsizes = quantile_bins(costs, K)
        N_samp = sampled_homogeneous_distribution(Float64.(labels), K - 1, n;
                     samples_per_class=S_N, rng=MersenneTwister(seed + 777))

        tv = true_vals(costs, n, scheds)
        gridceil = maximum(tv) / c_opt

        fP_ramp(θ) = binned_proxy_value(N_samp, binmeans, binsizes, n,
                                        linear_ramp(θ[1], θ[2], θ[3], θ[4], P)...)
        inits = Dict("transfer_1" => transfer_idx[1], "universal" => universal_idx)
        for (init, idx0) in inits
            c0 = combos[idx0]
            push!(rows, join(Any[fam, n, inst, seed, m, law, init, "init",
                                 tv[idx0] / c_opt, gridceil - tv[idx0] / c_opt], ","))
            θr, _ = compass_search(fP_ramp, [c0...],
                                   [(RAMP_γ̂[2]-RAMP_γ̂[1])/w̄/2, (RAMP_γ̂[2]-RAMP_γ̂[1])/w̄/2,
                                    (RAMP_β[2]-RAMP_β[1])/2, (RAMP_β[2]-RAMP_β[1])/2])
            γr, βr = linear_ramp(θr[1], θr[2], θr[3], θr[4], P)
            ar = qaoa_expectation(costs, n, γr, βr) / c_opt
            push!(rows, join(Any[fam, n, inst, seed, m, law, init, "polished",
                                 ar, gridceil - ar], ","))
            println("  inst=$inst $init: init=$(round(tv[idx0]/c_opt; digits=4)) polished=$(round(ar; digits=4))")
        end
        flush(stdout)
    end

    suffix = SMOKE ? "_smoke" : "_task$tid"
    outpath = joinpath(@__DIR__, "results$suffix.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,law,init,stage,ar,regret_vs_gridceil")
        foreach(r -> println(io, r), rows)
    end
    println("E031 task complete → $outpath")
end

main()
