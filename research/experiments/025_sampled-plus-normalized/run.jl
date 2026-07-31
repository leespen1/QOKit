#=
E025 — Does the recommended recipe compose? Sampled N + normalized objective.

The paper's decision rule recommends, in the statevector-affordable regime,
"sampled N (S ~ 10) with the empirical cost distribution and the normalized
objective at depth." But its two ingredients were validated separately:
E010 tested sampled N under the UNNORMALIZED convention, and E020/E024
tested normalization with the EXACT N. This experiment tests the actual
composition on the standard instance set (7 families x n in {12,14} x 10,
p = 1 and 3): regret under the four combinations
{exact, sampled-S10} x {unnormalized, normalized}.

Run: JULIA_NUM_THREADS=auto julia --project research/experiments/025_sampled-plus-normalized/run.jl (~15 min)
Smoke: E25_SMOKE=1 ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E25_SMOKE", "0") == "1"
const SEED  = 20260611
const NS    = SMOKE ? [10] : [12, 14]
const INST  = SMOKE ? 2 : 10
const GLEN  = SMOKE ? 10 : 40
const RLEN  = SMOKE ? 4 : 8
const S_N   = 10

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

"Both objectives (unnormalized, normalized) for every schedule column."
function objectives(N, γmat, βmat, counts, copt, n)
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
    jobs = [(fi, fam, n, inst) for (fi, (fam, _)) in enumerate(FAMILIES)
            for n in NS for inst in 1:INST]
    prep = map(jobs) do (fi, fam, n, inst)
        seed = instance_seed(fi, n, inst)
        edges = FAMILIES[fi][2](MersenneTwister(seed), n)
        m = length(edges); costs = maxcut_costs(n, edges); copt = maximum(costs)
        Nex = get_homogeneous_distribution_from_costs_direct(costs, m, n)
        Nsamp = sampled_homogeneous_distribution(
            costs, m, n; samples_per_class=S_N, rng=MersenneTwister(seed + 777))
        counts = zeros(Int, m + 1); for c in costs; counts[Int(c)+1] += 1; end
        (; fam, n, inst, m, costs, copt, Nex, Nsamp, counts)
    end
    println("prep done ($(length(prep)) instances)")
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
            regs = Float64[]
            for N in (q.Nex, q.Nsamp)
                fraw, fnorm = objectives(N, γmat, βmat, q.counts, q.copt, q.n)
                push!(regs, ceilv - F[argmax(fraw)])
                push!(regs, ceilv - F[argmax(fnorm)])
            end
            push!(out, join((q.fam, q.n, q.inst, p, round(ceilv, digits=6),
                round(regs[1], digits=6), round(regs[2], digits=6),
                round(regs[3], digits=6), round(regs[4], digits=6)), ","))
        end
        rows[k] = out
        println("done: ", q.fam, " n=", q.n, " inst=", q.inst)
    end
    path = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,p,ceiling," *
                    "regret_exact_raw,regret_exact_norm,regret_samp_raw,regret_samp_norm")
        for rs in rows, r in rs; println(io, r); end
    end
    println("wrote $path")
end

main()
