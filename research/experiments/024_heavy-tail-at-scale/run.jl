#=
E024 — Does the heavy-tail crossover survive scale?

E022/E023 established the transfer↔proxy inversion at the infinite-variance
boundary (α ≈ 2) at n=14–16. Here: n ∈ {22, 26}, α ∈ {1.5, 3.0} (one rung
each side of the boundary), p=1. At these sizes the exact binned homodist
is O(4ⁿ)-impossible, so the proxy side is the SAMPLED binned variant
(K=64, S=10 per bin) — the practical method whose niche is being claimed.
Ceilings by GPU statevector sweep on each instance's mean-weight-scaled
grid; transfer = grid-index transfer from weighted n=10 sources (the
practitioner's rescaling), as in E022/E023.

Tasks 1–8 = (ER(0.5), 3-regular) × (n = 22, 26) × (α = 1.5, 3.0),
5 instances each.
Submit:  cd research/experiments/024_heavy-tail-at-scale && sbatch run.sb
Smoke:   E24_SMOKE=1 julia --project research/experiments/024_heavy-tail-at-scale/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister, rand
using Base.Threads: @threads
using Statistics: mean

const USE_GPU = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end
println("USE_GPU = ", USE_GPU)

const SMOKE = get(ENV, "E24_SMOKE", "0") == "1"
SMOKE || USE_GPU || error("E024 requires a functional GPU (smoke runs on CPU)")

const SEED = 20260712
const SRC_SEED = 20260713
const SOURCE_N = 10
const N_SOURCES = SMOKE ? 1 : 3
const K = 64
const S_N = 10
const NINST = SMOKE ? 2 : 5
const GRID_LEN = SMOKE ? 10 : 40
const P1_γ̂ = collect(range(0.0, 2π; length=GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=GRID_LEN))

const TASKS = [(f, n, α) for α in (1.5, 3.0) for n in (22, 26),
                             f in ("ER(0.5)", "3-regular")] |> vec
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

"Allocation-free binned proxy sweep; returns argmax index over schedules."
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

"True ⟨C⟩ at one p=1 schedule (GPU when available)."
function true_exp(costs, costs_dev, n, γs, βs)
    if USE_GPU
        return gpu_qaoa_expectation_batched(costs_dev, n, γs, βs)
    end
    return qaoa_expectation(costs, n, γs, βs)
end

"CPU source landscape values (n=10)."
function src_vals(costs, scheds)
    vals = zeros(length(scheds))
    @threads for k in eachindex(scheds)
        γs, βs = scheds[k]
        vals[k] = qaoa_expectation(costs, SOURCE_N, γs, βs)
    end
    vals
end

function main()
    tid = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
    fam, n, α = SMOKE ? ("ER(0.5)", 12, 1.5) : TASKS[tid == 0 ? 1 : tid]
    draw = r -> (1 - r)^(-1 / α)
    println("task $tid: $fam n=$n alpha=$α K=$K S=$S_N")

    scheds_for(w̄) = vec([([γ̂ / w̄], [β]) for γ̂ in P1_γ̂, β in P1_β])

    # weighted sources on their own scaled grids (index ≡ shared γ̂ point)
    srcv = Dict{String, Vector{Vector{Float64}}}()
    for f in keys(GENS)
        reps = Vector{Vector{Float64}}()
        for rep in 1:N_SOURCES
            srng = MersenneTwister(SRC_SEED + 1000 * (f == "ER(0.5)" ? 1 : 2) + rep)
            sedges = GENS[f](srng, SOURCE_N)
            sw = draw.(rand(srng, length(sedges)))
            scosts = weighted_costs(SOURCE_N, sedges, sw)
            push!(reps, src_vals(scosts, scheds_for(mean(sw))) ./ maximum(scosts))
        end
        srcv[f] = reps
    end
    transfer_idx = [argmax(srcv[fam][r]) for r in 1:N_SOURCES]
    universal_idx = argmax(reduce(+, vcat([srcv[f] for f in keys(GENS)]...)))

    rows = String[]
    for inst in 1:NINST
        seed = SEED + 100 * n + 10_000 * round(Int, 10α) + inst
        rng = MersenneTwister(seed)
        edges = GENS[fam](rng, n)
        m = length(edges)
        w = draw.(rand(rng, m))
        costs = weighted_costs(n, edges, w)
        c_opt = maximum(costs)
        scheds = scheds_for(mean(w))
        costs_dev = USE_GPU ? CUDA.CuArray(costs) : nothing

        println("inst=$inst m=$m: bins...")
        flush(stdout)
        labels, binmeans, binsizes = quantile_bins(costs, K)
        N_samp = sampled_homogeneous_distribution(Float64.(labels), K - 1, n;
                     samples_per_class=S_N, rng=MersenneTwister(seed + 777))

        println("inst=$inst: ceiling...")
        flush(stdout)
        best = -Inf
        for (k, (γs, βs)) in enumerate(scheds)
            best = max(best, true_exp(costs, costs_dev, n, γs, βs))
            USE_GPU && k % 256 == 0 && CUDA.reclaim()
        end
        ceil_ar = best / c_opt

        emit(method, idx) = begin
            γs, βs = scheds[idx]
            ar = true_exp(costs, costs_dev, n, γs, βs) / c_opt
            @assert ceil_ar - ar > -1e-8
            push!(rows, join(Any[fam, n, inst, seed, m, α, method, ar,
                                 ceil_ar, ceil_ar - ar], ","))
        end
        push!(rows, join(Any[fam, n, inst, seed, m, α, "ceiling", ceil_ar,
                             ceil_ar, 0.0], ","))
        emit("binned_sampled", binned_proxy_argmax(N_samp, binmeans, binsizes, n, scheds))
        for r in 1:N_SOURCES
            emit("transfer_$r", transfer_idx[r])
        end
        emit("universal", universal_idx)
        println("done: inst=$inst ceil=$(round(ceil_ar; digits=4))")
        flush(stdout)
        USE_GPU && (costs_dev = nothing; CUDA.reclaim())
    end

    suffix = SMOKE ? "_smoke" : "_task$tid"
    outpath = joinpath(@__DIR__, "results$suffix.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,alpha,method,ar,ceil,regret")
        foreach(r -> println(io, r), rows)
    end
    println("E024 task complete → $outpath")
end

main()
