#=
E026 — The dispersion theory's prediction, tested on Pareto.

E025 unified the crossover as order-unity weight dispersion (CV²).
On Pareto(α), CV² = 1/(α(α−2)): 2.27 at α=2.2, 1.04 at α=2.4, 0.64 at
α=2.6. PRE-REGISTERED PREDICTION (before this run): the transfer↔proxy
flip sits between α=2.2 (proxy side) and α=2.4–2.6 (transfer side) —
i.e. inside the gap E023's grid skipped. Same harness as E023.

Output: results_task<ID>.csv (family,n,instance,seed,m,law,method,p,ar,ceil,regret).
Submit:  cd research/experiments/026_predicted-flip && sbatch run.sb
Smoke:   E26_SMOKE=1 julia --project research/experiments/026_predicted-flip/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister, rand
using Base.Threads: @threads
using Statistics: mean

const SMOKE = get(ENV, "E26_SMOKE", "0") == "1"

const SEED = 20260716
const SRC_SEED = 20260717
const SOURCE_N = 10
const N_SOURCES = SMOKE ? 1 : 3
const K = 64
const S_N = 10
const NINST = SMOKE ? 2 : 10
const GRID_LEN = SMOKE ? 10 : 40
const RAMP_LEN = SMOKE ? 4 : 8
const P_RAMP = 3

const P1_γ̂ = collect(range(0.0, 2π; length=GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=GRID_LEN))
const RAMP_γ̂ = collect(range(0.1, 3.2; length=RAMP_LEN))
const RAMP_β = collect(range(0.05, 0.8; length=RAMP_LEN))

const ALPHAS = Dict("a2.2" => 2.2, "a2.4" => 2.4, "a2.6" => 2.6)
const LAWS = Dict(k => (r -> (1 - r)^(-1 / α)) for (k, α) in ALPHAS)
const TASKS = [(f, 16, l) for l in sort(collect(keys(ALPHAS)))
                           for f in ("ER(0.5)", "3-regular")] |> vec
const GENS = Dict(
    "ER(0.5)"   => (rng, n) -> erdos_renyi_edges(n, 0.5; rng),
    "3-regular" => (rng, n) -> random_regular_edges(n, 3; rng),
)

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
    sums = zeros(K)
    sizes = zeros(Int, K)
    for i in 1:N
        b = labels[i] + 1
        sums[b] += costs[i]
        sizes[b] += 1
    end
    @assert all(>(0), sizes)
    labels, sums ./ sizes, sizes
end

"Allocation-free binned proxy sweep (E020's kernel); returns argmax index."
function binned_proxy_argmax(N, binmeans, binsizes, n, schedules)
    K = length(binmeans)
    βset = Dict{Float64, Vector{ComplexF64}}()
    for (γs, βs) in schedules, β in βs
        haskey(βset, β) || (βset[β] = vec(get_β_factors([β], n)))
    end
    vals = zeros(length(schedules))
    nt = Threads.nthreads()
    @threads :static for t in 1:nt
        Q = zeros(ComplexF64, K)
        Qsrc = zeros(ComplexF64, K)
        Qnew = zeros(ComplexF64, K)
        for k in t:nt:length(schedules)
            γs, βs = schedules[k]
            fill!(Q, ComplexF64(1 / sqrt(2.0^n)))
            for ℓ in eachindex(γs)
                f = βset[βs[ℓ]]
                γ = γs[ℓ]
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

"True-landscape values over schedules (threaded); returns the vals vector."
function true_vals(costs, n, schedules)
    vals = zeros(length(schedules))
    @threads for k in eachindex(schedules)
        γs, βs = schedules[k]
        vals[k] = qaoa_expectation(costs, n, γs, βs)
    end
    vals
end

function main()
    tid = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
    fam, n, law = SMOKE ? ("ER(0.5)", 12, "a2.4") : TASKS[tid == 0 ? 1 : tid]
    draw = LAWS[law]
    println("task $tid: $fam n=$n law=$law K=$K S=$S_N")

    "Schedules on this instance's scale: γ = γ̂ / w̄."
    function scheds_for(w̄)
        p1 = vec([([γ̂ / w̄], [β]) for γ̂ in P1_γ̂, β in P1_β])
        combos = vec([(g1 / w̄, gf / w̄, b1, bf) for g1 in RAMP_γ̂, gf in RAMP_γ̂,
                                                    b1 in RAMP_β, bf in RAMP_β])
        p1, [linear_ramp(c..., P_RAMP) for c in combos]
    end

    # weighted sources on their OWN scaled grids: index k ≡ shared γ̂ point
    src_vals = Dict{String, Vector{Dict{Int, Vector{Float64}}}}()
    for f in keys(GENS)
        reps = Vector{Dict{Int, Vector{Float64}}}()
        for rep in 1:N_SOURCES
            srng = MersenneTwister(SRC_SEED + 1000 * (f == "ER(0.5)" ? 1 : 2) + rep)
            sedges = GENS[f](srng, SOURCE_N)
            sw = draw.(rand(srng, length(sedges)))
            scosts = weighted_costs(SOURCE_N, sedges, sw)
            c_opt = maximum(scosts)
            sp1, sp3 = scheds_for(mean(sw))
            push!(reps, Dict(1 => true_vals(scosts, SOURCE_N, sp1) ./ c_opt,
                             3 => true_vals(scosts, SOURCE_N, sp3) ./ c_opt))
        end
        src_vals[f] = reps
    end
    transfer_idx = Dict(p => [argmax(src_vals[fam][r][p]) for r in 1:N_SOURCES]
                        for p in (1, 3))
    universal_idx = Dict(p => argmax(reduce(+, [reps[p] for f in keys(GENS)
                                                 for reps in src_vals[f]]))
                         for p in (1, 3))

    rows = String[]
    for inst in 1:NINST
        seed = SEED + 100 * n + inst
        rng = MersenneTwister(seed)
        edges = GENS[fam](rng, n)
        m = length(edges)
        w = draw.(rand(rng, m))
        costs = weighted_costs(n, edges, w)
        p1_scheds, ramp_scheds = scheds_for(mean(w))
        c_opt = maximum(costs)
        labels, binmeans, binsizes = quantile_bins(costs, K)
        N_exact = get_homogeneous_distribution_from_costs_direct(Float64.(labels), K - 1, n)
        N_samp = sampled_homogeneous_distribution(Float64.(labels), K - 1, n;
                     samples_per_class=S_N, rng=MersenneTwister(seed + 777))

        for (p, scheds) in ((1, p1_scheds), (P_RAMP, ramp_scheds))
            tv = true_vals(costs, n, scheds)
            ceil_ar = maximum(tv) / c_opt
            emit(method, idx) = begin
                ar = tv[idx] / c_opt
                @assert ceil_ar - ar > -1e-8
                push!(rows, join(Any[fam, n, inst, seed, m, law, method, p, ar,
                                     ceil_ar, ceil_ar - ar], ","))
            end
            push!(rows, join(Any[fam, n, inst, seed, m, law, "ceiling", p, ceil_ar,
                                 ceil_ar, 0.0], ","))
            emit("binned_exact", binned_proxy_argmax(N_exact, binmeans, binsizes, n, scheds))
            emit("binned_sampled", binned_proxy_argmax(N_samp, binmeans, binsizes, n, scheds))
            for r in 1:N_SOURCES
                emit("transfer_$r", transfer_idx[p][r])
            end
            emit("universal", universal_idx[p])
        end
        println("done: inst=$inst m=$m")
        flush(stdout)
    end

    suffix = SMOKE ? "_smoke" : "_task$tid"
    outpath = joinpath(@__DIR__, "results$suffix.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,law,method,p,ar,ceil,regret")
        foreach(r -> println(io, r), rows)
    end
    println("E026 task complete → $outpath")
end

main()
