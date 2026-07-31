#=
E046 — Can a polynomial-time Monte Carlo estimator of N(c'; d, c) replace the
exponential-cost enumeration sampler for QAOA parameter setting?

Question: The paper's sampled-N recipe (E041/E044: sampled N, S=10, plus the
normalized objective) still enumerates all 2^n bitstrings, capping it at
n <= 30. Here we estimate everything with polynomial cost in n: Wang-Landau
for the density of states / class counts, a flat-histogram chain for class
members, and C(n,d)-scaled random d-subset sampling for the profile rows
(estimator.jl). Does {poly-estimated N + estimated counts + normalized
objective} match {enumeration-sampled S=10 + normalized} and
{exact + normalized} in regret (within ~0.01 AR)?

Design: 7 families x n in {12,14} x 5 instances (E041 seed scheme);
variants {exact, samp10, poly} x {raw, norm}; p=1 grid (40x40) and p=3
linear-ramp grid (8^4); regret vs the true statevector grid ceiling. Each
variant pairs N with its own counts (exact counts for exact/samp10,
Wang-Landau counts for poly). Also: wall-clock of poly vs enumeration
estimator per instance, plus a scaling timing run at n=20 (3 ER(0.5)
instances, estimator only).

Run:   JULIA_NUM_THREADS=auto julia --project research/experiments/046_polytime-N-estimator/run.jl  (~30-45 min)
Smoke: E46_SMOKE=1 julia --project research/experiments/046_polytime-N-estimator/run.jl  (~2 min)

Outputs: results.csv (regrets), timing.csv (estimator wall-clocks).
=#

using JuliaQAOA
using Random: MersenneTwister, Xoshiro
using Base.Threads: @threads

include(joinpath(@__DIR__, "estimator.jl"))

const SMOKE = get(ENV, "E46_SMOKE", "0") == "1"
const SEED  = 20260611
const NVALS = SMOKE ? [10] : [12, 14]
const INST  = SMOKE ? 2 : 5
const GLEN  = SMOKE ? 10 : 40
const RLEN  = SMOKE ? 4 : 8
const KSUB  = SMOKE ? 60 : 200
const WLKW  = SMOKE ? (; lnf_final=1e-2, stage_cap=200_000) : (; lnf_final=1e-3, stage_cap=2_000_000)
const SCALE_N = 20            # scaling timing run (full mode only)
const SCALE_INST = 3

const FAMILIES = [
    ("ER(0.5)",       (rng, n) -> erdos_renyi_edges(n, 0.5; rng)),
    ("ER(0.25)",      (rng, n) -> erdos_renyi_edges(n, 0.25; rng)),
    ("BA(k=2)",       (rng, n) -> barabasi_albert_edges(n, 2; rng)),
    ("BA(k=4)",       (rng, n) -> barabasi_albert_edges(n, 4; rng)),
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("WS(k=4;b=0.5)", (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]
instance_seed(fi, n, inst) = SEED + 10_000 * fi + 100 * n + inst

const P1γ = collect(range(0.0, π;   length = GLEN))
const P1β = collect(range(0.0, π/2; length = GLEN))
const Rγ  = collect(range(0.05, 1.6; length = RLEN))
const Rβ  = collect(range(0.05, 0.8; length = RLEN))
const P1_sched = vec([([g], [b]) for g in P1γ, b in P1β])
const R_combos = vec([(g1,gf,b1,bf) for g1 in Rγ, gf in Rγ, b1 in Rβ, bf in Rβ])
const R_sched  = [linear_ramp(c..., 3) for c in R_combos]

function objectives(N, γmat, βmat, counts::AbstractVector, copt)
    Q = QAOA_proxy_multi(N, γmat, βmat)[end]
    K = size(Q, 2)
    raw = zeros(Float64, K); nrm = zeros(Float64, K)
    for k in 1:K, ci in axes(Q, 1)
        w = counts[ci] * abs2(Q[ci, k])
        raw[k] += w * (ci - 1); nrm[k] += w
    end
    return raw ./ copt, raw ./ nrm ./ copt
end

function prepare(fi, fam, n, inst)
    seed = instance_seed(fi, n, inst)
    edges = FAMILIES[fi][2](MersenneTwister(seed), n)
    m = length(edges)
    costs = maxcut_costs(n, edges)
    copt = maximum(costs)
    Nex = get_homogeneous_distribution_from_costs_direct(costs, m, n)
    counts = zeros(Float64, m + 1)
    for c in costs; counts[Int(c)+1] += 1.0; end
    t_enum = @elapsed s10 = sampled_homogeneous_distribution(costs, m, n;
        samples_per_class=10, rng=MersenneTwister(seed + 777))
    t_poly = @elapsed begin
        Npoly, cpoly, members = polytime_homogeneous_distribution(n, edges;
            samples_per_class=10, K_subsets=KSUB, rng=Xoshiro(seed + 999), WLKW...)
    end
    # Fail fast if the sampler claims a class no bitstring attains
    # (would indicate a delta-cost bug); merely missing tiny classes is a
    # recorded outcome, not an error.
    attained = findall(>(0.0), counts)
    visited  = findall(ci -> cpoly[ci] > 0.0, 1:(m + 1))
    @assert issubset(visited, attained) "poly estimator visited an unattainable cost class"
    missed = length(setdiff(attained, visited))
    (; fam, n, inst, m, costs, copt, missed, t_enum, t_poly,
       variants=[("exact",  Nex,   counts),
                 ("samp10", s10,   counts),
                 ("poly",   Npoly, cpoly)])
end

function main()
    jobs = [(fi, fam, n, inst) for n in NVALS
            for (fi, (fam, _)) in enumerate(FAMILIES) for inst in 1:INST]
    prep = map(jobs) do (fi, fam, n, inst)
        q = prepare(fi, fam, n, inst)
        println("prep: $fam n=$n inst=$inst m=$(q.m) missed=$(q.missed) ",
                "t_enum=$(round(q.t_enum, digits=2))s t_poly=$(round(q.t_poly, digits=2))s")
        q
    end
    rows = Vector{Vector{String}}(undef, length(prep))
    @threads for k in eachindex(prep)
        q = prep[k]
        out = String[]
        for (p, scheds, γmat, βmat) in (
                (1, P1_sched,
                    reshape([s[1][1] for s in P1_sched], :, 1),
                    reshape([s[2][1] for s in P1_sched], :, 1)),
                (3, R_sched,
                    linear_ramp_matrix([c[1] for c in R_combos], [c[2] for c in R_combos],
                                       [c[3] for c in R_combos], [c[4] for c in R_combos], 3)...))
            F = [qaoa_expectation(q.costs, q.n, γs, βs) / q.copt for (γs, βs) in scheds]
            ceilv = maximum(F)
            for (kind, N, counts) in q.variants
                fraw, fnorm = objectives(N, γmat, βmat, counts, q.copt)
                push!(out, join((q.fam, q.n, q.inst, p, round(ceilv, digits=6), kind,
                    round(ceilv - F[argmax(fraw)], digits=6),
                    round(ceilv - F[argmax(fnorm)], digits=6)), ","))
            end
        end
        rows[k] = out
        println("done: ", q.fam, " n=", q.n, " inst=", q.inst)
    end
    respath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(respath, "w") do io
        println(io, "family,n,inst,p,ceiling,kind,regret_raw,regret_norm")
        for rs in rows, r in rs; println(io, r); end
    end
    println("wrote $respath")

    timpath = joinpath(@__DIR__, SMOKE ? "timing_smoke.csv" : "timing.csv")
    open(timpath, "w") do io
        println(io, "family,n,inst,m,t_enum_s,t_poly_s,classes_missed")
        for q in prep
            println(io, join((q.fam, q.n, q.inst, q.m,
                round(q.t_enum, digits=3), round(q.t_poly, digits=3), q.missed), ","))
        end
        if !SMOKE
            for inst in 1:SCALE_INST
                seed = instance_seed(1, SCALE_N, inst)
                edges = erdos_renyi_edges(SCALE_N, 0.5; rng=MersenneTwister(seed))
                m = length(edges)
                t_enum = @elapsed begin
                    costs = maxcut_costs(SCALE_N, edges)
                    sampled_homogeneous_distribution(costs, m, SCALE_N;
                        samples_per_class=10, rng=MersenneTwister(seed + 777))
                end
                t_poly = @elapsed begin
                    _, cpoly, _ = polytime_homogeneous_distribution(SCALE_N, edges;
                        samples_per_class=10, K_subsets=KSUB, rng=Xoshiro(seed + 999), WLKW...)
                end
                println(io, join(("ER(0.5)-scale", SCALE_N, inst, m,
                    round(t_enum, digits=3), round(t_poly, digits=3), ""), ","))
                println("scale: n=$SCALE_N inst=$inst m=$m t_enum=$(round(t_enum, digits=2))s ",
                        "t_poly=$(round(t_poly, digits=2))s")
            end
        end
    end
    println("wrote $timpath")
end

main()
