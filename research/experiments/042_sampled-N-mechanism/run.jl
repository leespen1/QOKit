#=
E042 — Why does sampled N beat exact N? Noise vs estimator bias.

E041 found that sampled N (S=10) + normalized objective beats exact N +
normalized at both depths. The sampler averages exact n(x;d,c) rows over
up to S members per cost class, so it is an unbiased estimator of the class
mean; any systematic gain must come either from the *noise* itself
(interacting with argmax selection) or from the small-class regime where
S >= class size makes the estimate exact. Three probes on n=12:

  A. R=10 independent sampling replicates per instance: the spread of
     per-replicate regret, and the regret of the replicate-AVERAGED N
     (noise suppressed, any bias kept). If averaged-N regresses to the
     exact-N regret, the active ingredient is noise.
  B. S-sweep {3, 10, 30, 100, 300}: if noise is the ingredient, the gain
     should fade as S grows and the sampled N converges to exact.

Run: JULIA_NUM_THREADS=auto julia --project research/experiments/042_sampled-N-mechanism/run.jl (~20 min)
Smoke: E26_SMOKE=1 ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E26_SMOKE", "0") == "1"
const SEED  = 20260611             # same instance seeds as E041
const N_V   = 12
const INST  = SMOKE ? 2 : 10
const GLEN  = SMOKE ? 10 : 40
const RLEN  = SMOKE ? 4 : 8
const REPS  = SMOKE ? 3 : 10
const SVALS = SMOKE ? [3, 30] : [3, 10, 30, 100, 300]

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
        # variants: (kind, param, N). Replicate rng streams are disjoint by
        # construction: rep r uses seed+1000r+777, sweep S uses seed+13S+778.
        variants = Vector{Tuple{String,Int,Array{Float64,3}}}()
        push!(variants, ("exact", 0, q.Nex))
        reps = [sampled_homogeneous_distribution(q.costs, q.m, N_V;
                    samples_per_class=10, rng=MersenneTwister(q.seed + 1000r + 777))
                for r in 1:REPS]
        for (r, Ns) in enumerate(reps)
            push!(variants, ("rep", r, Ns))
        end
        push!(variants, ("avg", REPS, sum(reps) ./ REPS))
        for S in SVALS
            push!(variants, ("sweep", S, sampled_homogeneous_distribution(
                q.costs, q.m, N_V; samples_per_class=S,
                rng=MersenneTwister(q.seed + 13S + 778))))
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
            for (kind, param, N) in variants
                fraw, fnorm = objectives(N, γmat, βmat, q.counts, q.copt)
                push!(out, join((q.fam, N_V, q.inst, p, round(ceilv, digits=6),
                    kind, param,
                    round(ceilv - F[argmax(fraw)], digits=6),
                    round(ceilv - F[argmax(fnorm)], digits=6)), ","))
            end
        end
        rows[k] = out
        println("done: ", q.fam, " inst=", q.inst)
    end
    path = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,p,ceiling,kind,param,regret_raw,regret_norm")
        for rs in rows, r in rs; println(io, r); end
    end
    println("wrote $path")
end

main()
