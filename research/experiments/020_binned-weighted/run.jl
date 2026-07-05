#=
E020 — Binned compression for WEIGHTED MaxCut: does the calculus survive
continuous costs?

With generic real edge weights every cost class is a singleton, the
projector is the identity, and the integer-cost compression trivializes —
the paper's top referee risk. But Theorems 2–3 hold for ANY fixed partition
of the computational basis. Here the partition is K equal-population
(quantile) bins of the weighted cost, and we measure:

  (i)  one-layer leakage λ(γ,β; K) from |+⟩ — expected to decompose into
       the structural compression term plus a binning term that shrinks as
       the within-bin cost spread ~ 1/K;
  (ii) parameter-setting regret of the binned proxy (bin-label homodist via
       the existing integer machinery; bin-mean phases; expectation with
       bin-mean values) against the true weighted grid ceiling, at p=1 and
       p=3 ramps, for K ∈ {4, 8, 16, 32, 64, 128}.

Weights are i.i.d. U[0,1] (fixed seeds). Because mean weight is 1/2, the
natural γ scale roughly doubles vs the unweighted grid: p=1 uses
γ ∈ [0, 2π]; ramp endpoints γ ∈ [0.1, 3.2].

Slurm array: tasks 1–4 = (ER(0.5), 3-regular) × (n = 14, 16), 10 instances.
Submit:  cd research/experiments/020_binned-weighted && sbatch run.sb
Smoke:   E20_SMOKE=1 julia --project research/experiments/020_binned-weighted/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister, rand
using Base.Threads: @threads
using Statistics: mean

const SMOKE = get(ENV, "E20_SMOKE", "0") == "1"

const SEED = 20260706
const KS = SMOKE ? [4, 16] : [4, 8, 16, 32, 64, 128]
const NINST = SMOKE ? 2 : 10
const GRID_LEN = SMOKE ? 10 : 40
const RAMP_LEN = SMOKE ? 4 : 8
const P_RAMP = 3
const LEAK_ANGLES = [(0.2, 0.2), (0.5, 0.3), (1.0, 0.4)]

const P1_γ = collect(range(0.0, 2π; length=GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=GRID_LEN))
const RAMP_γ = collect(range(0.1, 3.2; length=RAMP_LEN))
const RAMP_β = collect(range(0.05, 0.8; length=RAMP_LEN))

const TASKS = [("ER(0.5)", 14), ("3-regular", 14), ("ER(0.5)", 16), ("3-regular", 16)]
const GENS = Dict(
    "ER(0.5)"   => (rng, n) -> erdos_renyi_edges(n, 0.5; rng),
    "3-regular" => (rng, n) -> random_regular_edges(n, 3; rng),
)

"Weighted cut costs, threaded."
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

"Equal-population bins: 0-based labels, bin means, bin sizes."
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
    @assert all(>(0), sizes) "empty bin (K too large for distinct values?)"
    labels, sums ./ sizes, sizes
end

"One-layer leakage from |+⟩ under the TRUE weighted layer, onto the bins."
function one_layer_leakage(costs, labels, K, n, γ, β)
    ψ = fill(ComplexF64(1 / sqrt(2.0^n)), 1 << n)
    apply_phase_gate!(ψ, costs, γ)
    apply_x_mixer!(ψ, β, n)
    sums = zeros(ComplexF64, K)
    sizes = zeros(Int, K)
    for i in eachindex(ψ)
        b = labels[i] + 1
        sums[b] += ψ[i]
        sizes[b] += 1
    end
    means = sums ./ sizes
    resid2 = 0.0
    for i in eachindex(ψ)
        resid2 += abs2(ψ[i] - means[labels[i] + 1])
    end
    sqrt(resid2)
end

#=
Binned proxy sweep. State Q ∈ C^K is the shared per-bitstring amplitude of
each bin (paper convention). One layer:
  Q'(b') = Σ_{d,b} f_d(β) · cis(−γ·c̄_b/2) · Q(b) · N[b', d, b]
Expectation (bin-mean values): ⟨C⟩ = Σ_b M_b |Q(b)|² c̄_b.
Returns the argmax index over the schedule list.
=#
function binned_proxy_argmax(N, binmeans, binsizes, n, schedules)
    K = length(binmeans)
    vals = zeros(length(schedules))
    @threads for k in eachindex(schedules)
        γs, βs = schedules[k]
        Q = fill(ComplexF64(1 / sqrt(2.0^n)), K)
        for ℓ in eachindex(γs)
            f = vec(get_β_factors([βs[ℓ]], n))          # (n+1)
            phase = [cis(-γs[ℓ] * c / 2) for c in binmeans]
            Qsrc = phase .* Q
            Qnew = zeros(ComplexF64, K)
            for b in 1:K, d in 0:n
                w = f[d + 1] * Qsrc[b]
                w == 0 && continue
                for b′ in 1:K
                    Qnew[b′] += w * N[b′, d + 1, b]
                end
            end
            Q = Qnew
        end
        vals[k] = sum(binsizes[b] * abs2(Q[b]) * binmeans[b] for b in 1:K)
    end
    argmax(vals)
end

"True ⟨C⟩ over schedules, threaded CPU; returns (ceiling, best index)."
function true_grid(costs, n, schedules)
    vals = zeros(length(schedules))
    @threads for k in eachindex(schedules)
        γs, βs = schedules[k]
        vals[k] = qaoa_expectation(costs, n, γs, βs)
    end
    maximum(vals)
end

function main()
    tid = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
    fam, n = SMOKE ? ("ER(0.5)", 12) : TASKS[tid == 0 ? 1 : tid]
    println("task $tid: $fam n=$n instances=$NINST Ks=$KS")

    p1_scheds = vec([([γ], [β]) for γ in P1_γ, β in P1_β])
    ramp_combos = vec([(g1, gf, b1, bf) for g1 in RAMP_γ, gf in RAMP_γ,
                                            b1 in RAMP_β, bf in RAMP_β])
    ramp_scheds = [linear_ramp(c..., P_RAMP) for c in ramp_combos]

    rows = String[]
    for inst in 1:NINST
        seed = SEED + 100 * n + inst
        rng = MersenneTwister(seed)
        edges = GENS[fam](rng, n)
        m = length(edges)
        w = rand(rng, m)
        costs = weighted_costs(n, edges, w)
        c_opt = maximum(costs)

        ceil1 = true_grid(costs, n, p1_scheds) / c_opt
        ceil3 = true_grid(costs, n, ramp_scheds) / c_opt
        push!(rows, join(Any[fam, n, inst, seed, m, "ceiling", 0, 1, ceil1, ceil1, 0.0], ","))
        push!(rows, join(Any[fam, n, inst, seed, m, "ceiling", 0, 3, ceil3, ceil3, 0.0], ","))

        for K in KS
            labels, binmeans, binsizes = quantile_bins(costs, K)
            for (γ, β) in LEAK_ANGLES
                λ = one_layer_leakage(costs, labels, K, n, γ, β)
                push!(rows, join(Any[fam, n, inst, seed, m, "leak_g$(γ)_b$(β)",
                                     K, 1, λ, "", ""], ","))
            end
            N = get_homogeneous_distribution_from_costs_direct(Float64.(labels), K - 1, n)
            for (p, scheds) in ((1, p1_scheds), (P_RAMP, ramp_scheds))
                best = binned_proxy_argmax(N, binmeans, binsizes, n, scheds)
                γs, βs = scheds[best]
                ar = qaoa_expectation(costs, n, γs, βs) / c_opt
                ceilp = p == 1 ? ceil1 : ceil3
                @assert ceilp - ar > -1e-8 "binned proxy beat the ceiling"
                push!(rows, join(Any[fam, n, inst, seed, m, "binned_proxy",
                                     K, p, ar, ceilp, ceilp - ar], ","))
            end
        end
        println("done: inst=$inst m=$m ceil1=$(round(ceil1; digits=4))")
        flush(stdout)
    end

    suffix = SMOKE ? "_smoke" : "_task$tid"
    outpath = joinpath(@__DIR__, "results$suffix.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,kind,K,p,value,ceil,regret")
        foreach(r -> println(io, r), rows)
    end
    println("E020 task complete → $outpath")
end

main()
