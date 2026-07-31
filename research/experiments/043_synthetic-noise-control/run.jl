#=
E043 — Synthetic-noise control: does generic noise on exact N reproduce the
sampled-N gain?

E042 showed the sampled-N advantage (under the normalized objective) is
carried by the sampling noise itself. Two hypotheses remain:
  (H1) generic stochastic regularization — any zero-mean perturbation of N
       of comparable size disrupts the winner's curse; then multiplicative
       iid noise on exact N should reproduce the gain, and subsampling is
       just one way to buy noise.
  (H2) structured heterogeneity — the sampler's noise directions are
       differences of REAL member profiles n(x;d,c), and that structure is
       essential; then iid noise should fail to help (or hurt).

Design: exact N perturbed entrywise, N'[i] = N[i] * max(0, 1 + sigma*z),
z ~ Normal(0,1) iid, sigma in {0.05, 0.1, 0.2, 0.4, 0.8}, R=5 replicates
per sigma; same 70 instances (7 families x n=12), both depths, both
objectives. Reference rows: exact N and one S=10 sampled draw (E042 stream).

Run: JULIA_NUM_THREADS=auto julia --project research/experiments/043_synthetic-noise-control/run.jl (~15 min)
Smoke: E27_SMOKE=1 ...
=#

using JuliaQAOA
using Random: MersenneTwister, randn
using Base.Threads: @threads

const SMOKE = get(ENV, "E27_SMOKE", "0") == "1"
const SEED  = 20260611
const N_V   = 12
const INST  = SMOKE ? 2 : 10
const GLEN  = SMOKE ? 10 : 40
const RLEN  = SMOKE ? 4 : 8
const REPS  = SMOKE ? 2 : 5
const SIGMAS = SMOKE ? [0.2] : [0.05, 0.1, 0.2, 0.4, 0.8]

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

function objectives(N, γmat, βmat, counts, copt)
    Q = QAOA_proxy_multi(N, γmat, βmat)[end]
    K = size(Q, 2)
    raw = zeros(Float64, K); nrm = zeros(Float64, K)
    for k in 1:K, ci in axes(Q, 1)
        w = counts[ci] * abs2(Q[ci, k])
        raw[k] += w * (ci - 1); nrm[k] += w
    end
    return raw ./ copt, raw ./ nrm ./ copt
end

perturb(N, σ, rng) = N .* max.(0.0, 1.0 .+ σ .* randn(rng, size(N)))

function main()
    jobs = [(fi, fam, inst) for (fi, (fam, _)) in enumerate(FAMILIES) for inst in 1:INST]
    prep = map(jobs) do (fi, fam, inst)
        seed = instance_seed(fi, N_V, inst)
        edges = FAMILIES[fi][2](MersenneTwister(seed), N_V)
        m = length(edges); costs = maxcut_costs(N_V, edges); copt = maximum(costs)
        Nex = get_homogeneous_distribution_from_costs_direct(costs, m, N_V)
        counts = zeros(Int, m + 1); for c in costs; counts[Int(c)+1] += 1; end
        (; fam, inst, seed, m, costs, copt, Nex, counts)
    end
    println("prep done ($(length(prep)) instances)")
    rows = Vector{Vector{String}}(undef, length(prep))
    @threads for k in eachindex(prep)
        q = prep[k]
        variants = Vector{Tuple{String,Float64,Int,Array{Float64,3}}}()
        push!(variants, ("exact", 0.0, 0, q.Nex))
        push!(variants, ("samp10", 0.0, 1, sampled_homogeneous_distribution(
            q.costs, q.m, N_V; samples_per_class=10,
            rng=MersenneTwister(q.seed + 1777))))
        for σ in SIGMAS, r in 1:REPS
            push!(variants, ("iid", σ, r, perturb(
                q.Nex, σ, MersenneTwister(q.seed + round(Int, 10_000σ) + 31r))))
        end
        out = String[]
        for (p, scheds, γmat, βmat) in (
                (1, P1_sched,
                    reshape([s[1][1] for s in P1_sched], :, 1),
                    reshape([s[2][1] for s in P1_sched], :, 1)),
                (3, R_sched,
                    linear_ramp_matrix([c[1] for c in R_combos], [c[2] for c in R_combos],
                                       [c[3] for c in R_combos], [c[4] for c in R_combos], 3)...))
            F = [qaoa_expectation(q.costs, N_V, γs, βs) / q.copt for (γs, βs) in scheds]
            ceilv = maximum(F)
            for (kind, σ, r, N) in variants
                fraw, fnorm = objectives(N, γmat, βmat, q.counts, q.copt)
                push!(out, join((q.fam, N_V, q.inst, p, round(ceilv, digits=6),
                    kind, σ, r,
                    round(ceilv - F[argmax(fraw)], digits=6),
                    round(ceilv - F[argmax(fnorm)], digits=6)), ","))
            end
        end
        rows[k] = out
        println("done: ", q.fam, " inst=", q.inst)
    end
    path = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,p,ceiling,kind,sigma,rep,regret_raw,regret_norm")
        for rs in rows, r in rs; println(io, r); end
    end
    println("wrote $path")
end

main()
