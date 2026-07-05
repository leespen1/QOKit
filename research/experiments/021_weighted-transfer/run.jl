#=
E021 — Does transfer's dominance survive random weights?

E018/E019 showed plain within-family angle transfer beating every proxy
variant on UNWEIGHTED families — where instances share exact parameter
concentration. Random edge weights individualize every instance: a source's
optimal angles need not transfer, while the binned proxy (E020) sees each
instance's own cost structure. This is the first setting where the proxy
could genuinely win; either outcome completes the practitioner map.

Cells: (ER(0.5), 3-regular) × n ∈ {14, 16}, 10 instances, i.i.d. U[0,1]
edge weights (E020's seed convention so ceilings are comparable).
Depths: p=1 (γ ∈ [0, 2π] × β ∈ [0, π/2], 40×40) and p=3 ramps
(γ ∈ [0.1, 3.2], β ∈ [0.05, 0.8], 8⁴) — E020's grids.

Methods per instance (same grid as the true weighted ceiling; regret ≥ 0):
  ceiling          — true weighted grid ceiling (threaded CPU statevector)
  binned_exact     — E020's binned proxy, K=64, exact bin-label homodist
  binned_sampled   — same, but N estimated from S=10 samples per bin (the
                     practical variant: O(S·K·2ⁿ) not O(4ⁿ))
  transfer_1..3    — angles grid-optimized on the true landscape of an
                     n=10 same-family instance WITH ITS OWN random weights
  universal        — pooled argmax over all 6 weighted sources

Output: results_task<ID>.csv (long; family,n,instance,seed,m,method,p,ar,ceil,regret).
Submit:  cd research/experiments/021_weighted-transfer && sbatch run.sb
Smoke:   E21_SMOKE=1 julia --project research/experiments/021_weighted-transfer/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister, rand
using Base.Threads: @threads
using Statistics: mean

const SMOKE = get(ENV, "E21_SMOKE", "0") == "1"

const SEED = 20260706            # E020's convention: same targets
const SRC_SEED = 20260707
const SOURCE_N = 10
const N_SOURCES = SMOKE ? 1 : 3
const K = 64
const S_N = 10
const NINST = SMOKE ? 2 : 10
const GRID_LEN = SMOKE ? 10 : 40
const RAMP_LEN = SMOKE ? 4 : 8
const P_RAMP = 3

const P1_γ = collect(range(0.0, 2π; length=GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=GRID_LEN))
const RAMP_γ = collect(range(0.1, 3.2; length=RAMP_LEN))
const RAMP_β = collect(range(0.05, 0.8; length=RAMP_LEN))

const TASKS = [("ER(0.5)", 14), ("3-regular", 14), ("ER(0.5)", 16), ("3-regular", 16)]
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
    fam, n = SMOKE ? ("ER(0.5)", 12) : TASKS[tid == 0 ? 1 : tid]
    println("task $tid: $fam n=$n K=$K S=$S_N")

    p1_scheds = vec([([γ], [β]) for γ in P1_γ, β in P1_β])
    ramp_combos = vec([(g1, gf, b1, bf) for g1 in RAMP_γ, gf in RAMP_γ,
                                            b1 in RAMP_β, bf in RAMP_β])
    ramp_scheds = [linear_ramp(c..., P_RAMP) for c in ramp_combos]

    # weighted sources: both families pooled for `universal`, own family for transfer
    src_vals = Dict{String, Vector{Dict{Int, Vector{Float64}}}}()
    for f in keys(GENS)
        reps = Vector{Dict{Int, Vector{Float64}}}()
        for rep in 1:N_SOURCES
            srng = MersenneTwister(SRC_SEED + 1000 * (f == "ER(0.5)" ? 1 : 2) + rep)
            sedges = GENS[f](srng, SOURCE_N)
            sw = rand(srng, length(sedges))
            scosts = weighted_costs(SOURCE_N, sedges, sw)
            c_opt = maximum(scosts)
            push!(reps, Dict(1 => true_vals(scosts, SOURCE_N, p1_scheds) ./ c_opt,
                             3 => true_vals(scosts, SOURCE_N, ramp_scheds) ./ c_opt))
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
        w = rand(rng, m)
        costs = weighted_costs(n, edges, w)
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
                push!(rows, join(Any[fam, n, inst, seed, m, method, p, ar,
                                     ceil_ar, ceil_ar - ar], ","))
            end
            push!(rows, join(Any[fam, n, inst, seed, m, "ceiling", p, ceil_ar,
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
        println(io, "family,n,instance,seed,m,method,p,ar,ceil,regret")
        foreach(r -> println(io, r), rows)
    end
    println("E021 task complete → $outpath")
end

main()
