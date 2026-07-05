#=
E029 — Is the deep high-dispersion corner a grid artifact or intrinsic?

E028 found that at (high weight dispersion × p=10–20 ramps) every method is
bad (regret 0.03–0.3 against the 8⁴ ramp-endpoint grid ceiling). Two
hypotheses: (a) the coarse endpoint grid is the bottleneck — continuous
refinement recovers most of the value for everyone; (b) the corner is
intrinsically hard — refinement from any method's choice stalls far from
the best known schedule.

Design: E028's two worst cells, same seeds/instances (5 each):
  task 1: pareto15, ER(0.5), n=16, p=20
  task 2: logn2.5, 3-regular, n=16, p=20
For each instance, take four grid choices — grid ceiling, binned_sampled
proxy, best transfer source, universal — and refine each on the TRUE
landscape two ways:
  ramp4   — compass search over the 4 ramp endpoints (continuous ramps);
  full2p  — compass search over all 2p angles, started from ramp4's result.
Reference V* = best value seen anywhere for that instance. Report each
(start, stage) as AR and regret vs V*.

Output: results_task<ID>.csv
(family,n,instance,seed,m,law,start,stage,ar,ar_vstar,regret_vstar)
Submit:  cd research/experiments/029_corner-refinement && sbatch run.sb
Smoke:   E29_SMOKE=1 julia --project research/experiments/029_corner-refinement/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister, rand, randn
using Base.Threads: @threads
using Statistics: mean

const SMOKE = get(ENV, "E29_SMOKE", "0") == "1"

const SEED = 20260720            # E028's convention: same instances
const SRC_SEED = 20260721
const SOURCE_N = 10
const N_SOURCES = SMOKE ? 1 : 3
const K = 64
const S_N = 10
const NINST = SMOKE ? 1 : 5
const RAMP_LEN = SMOKE ? 4 : 8
const P = SMOKE ? 5 : 20

const RAMP_γ̂ = collect(range(0.1, 3.2; length=RAMP_LEN))
const RAMP_β = collect(range(0.05, 0.8; length=RAMP_LEN))

const TASKS = [("ER(0.5)", "pareto15"), ("3-regular", "logn2.5")]
const N_TARGET = SMOKE ? 12 : 16
const GENS = Dict(
    "ER(0.5)"   => (rng, n) -> erdos_renyi_edges(n, 0.5; rng),
    "3-regular" => (rng, n) -> random_regular_edges(n, 3; rng),
)
const LAWS = Dict("pareto15" => (rng, m) -> (1 .- rand(rng, m)).^(-1/1.5),
                  "logn2.5"  => (rng, m) -> exp.(2.5 .* randn(rng, m)))

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
    fam, law = TASKS[tid == 0 ? 1 : tid]
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
        best_r = argmax([tv[transfer_idx[r]] for r in 1:N_SOURCES])
        starts = Dict(
            "grid_ceiling"   => argmax(tv),
            "binned_sampled" => binned_proxy_argmax(N_samp, binmeans, binsizes, n, scheds),
            "best_transfer"  => transfer_idx[best_r],
            "universal"      => universal_idx,
        )

        f_ramp(θ) = qaoa_expectation(costs, n, linear_ramp(θ[1], θ[2], θ[3], θ[4], P)...)
        f_full(θ) = qaoa_expectation(costs, n, θ[1:P], θ[(P + 1):(2P)])

        results = Dict{Tuple{String, String}, Float64}()
        vstar = maximum(tv)
        for (name, idx) in starts
            c = combos[idx]
            results[(name, "grid")] = tv[idx]
            θr, vr = compass_search(f_ramp, [c...],
                                    [step for step in ((RAMP_γ̂[2]-RAMP_γ̂[1])/w̄/2,
                                     (RAMP_γ̂[2]-RAMP_γ̂[1])/w̄/2,
                                     (RAMP_β[2]-RAMP_β[1])/2, (RAMP_β[2]-RAMP_β[1])/2)])
            results[(name, "ramp4")] = vr
            γs, βs = linear_ramp(θr[1], θr[2], θr[3], θr[4], P)
            θf, vf = compass_search(f_full, vcat(γs, βs),
                                    fill(0.05 / w̄, P) |> x -> vcat(x, fill(0.05, P)))
            results[(name, "full2p")] = vf
            vstar = max(vstar, vr, vf)
            println("  inst=$inst $name: grid=$(round(tv[idx]/c_opt; digits=4)) ramp4=$(round(vr/c_opt; digits=4)) full2p=$(round(vf/c_opt; digits=4))")
            flush(stdout)
        end
        for ((name, stage), v) in results
            push!(rows, join(Any[fam, n, inst, seed, m, law, name, stage,
                                 v / c_opt, vstar / c_opt, (vstar - v) / c_opt], ","))
        end
        println("done inst=$inst vstar=$(round(vstar/c_opt; digits=4))")
        flush(stdout)
    end

    suffix = SMOKE ? "_smoke" : "_task$tid"
    outpath = joinpath(@__DIR__, "results$suffix.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,law,start,stage,ar,ar_vstar,regret_vstar")
        foreach(r -> println(io, r), rows)
    end
    println("E029 task complete → $outpath")
end

main()
