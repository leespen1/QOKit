#=
E022 — Why is the two-point certificate 10x loose? The exact signed
decomposition of regret.

With the signed error field e(theta) = F(theta) - F-hat(theta) (true minus
normalized-proxy objective), regret decomposes EXACTLY:

  regret = F(theta*) - F(theta-hat)
         = [e(theta*) - e(theta-hat)]  -  [F-hat(theta-hat) - F-hat(theta*)]
         =        Delta-e              -            margin,

where margin >= 0 by definition of theta-hat. The E020 certificate replaces
Delta-e by |e(theta*)| + |e(theta-hat)| and drops the margin; if the signed
errors at the two points are large but similar (positively correlated across
instances and angles), Delta-e is tiny while the certificate's sum is not.
This experiment measures the decomposition on the standard 140 instances at
p = 1 and 3: the identity (asserted to 1e-9), the cancellation factor
Delta-e / (|e*| + |e-hat|), the across-instance correlation of e(theta*) and
e(theta-hat), and the share of regret absorbed by the margin.

Run: JULIA_NUM_THREADS=auto julia --project research/experiments/022_signed-regret-decomposition/run.jl  (~10 min)
Smoke: E22_SMOKE=1 ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E22_SMOKE", "0") == "1"
const SEED  = 20260611
const NS    = SMOKE ? [10] : [12, 14]
const INST  = SMOKE ? 2 : 10
const GLEN  = SMOKE ? 10 : 40
const RLEN  = SMOKE ? 4 : 8

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

function normalized_proxy_landscape(N, γmat, βmat, counts, copt)
    Qs = QAOA_proxy_multi(N, γmat, βmat)
    Q = Qs[end]
    K = size(Q, 2)
    out = zeros(Float64, K)
    for k in 1:K
        raw = 0.0; nrm = 0.0
        for ci in axes(Q, 1)
            w = counts[ci] * abs2(Q[ci, k])
            raw += w * (ci - 1); nrm += w
        end
        out[k] = raw / nrm / copt
    end
    return out
end

function instance_rows(fam, n, inst, m, costs, copt, N, counts)
    rows = String[]
    for (p, scheds, γmat, βmat) in (
            (1, P1_sched,
                reshape([s[1][1] for s in P1_sched], :, 1),
                reshape([s[2][1] for s in P1_sched], :, 1)),
            (3, R_sched,
                linear_ramp_matrix([c[1] for c in R_combos], [c[2] for c in R_combos],
                                   [c[3] for c in R_combos], [c[4] for c in R_combos], 3)...))
        F = [qaoa_expectation(costs, n, γs, βs) / copt for (γs, βs) in scheds]
        Fhat = normalized_proxy_landscape(N, γmat, βmat, counts, copt)
        istar = argmax(F); ihat = argmax(Fhat)
        regret = F[istar] - F[ihat]
        estar = F[istar] - Fhat[istar]
        ehat  = F[ihat]  - Fhat[ihat]
        de = estar - ehat
        margin = Fhat[ihat] - Fhat[istar]
        margin >= -1e-12 || error("negative margin: $fam $n $inst p=$p")
        abs(regret - (de - margin)) < 1e-9 ||
            error("decomposition identity fails: $fam $n $inst p=$p")
        cancel = de / max(abs(estar) + abs(ehat), 1e-12)
        push!(rows, join((fam, n, inst, m, p,
            round(regret, digits=6), round(estar, digits=6), round(ehat, digits=6),
            round(de, digits=6), round(margin, digits=6),
            round(cancel, digits=6)), ","))
    end
    return rows
end

function main()
    jobs = [(fi, fam, n, inst) for (fi, (fam, _)) in enumerate(FAMILIES)
            for n in NS for inst in 1:INST]
    prep = map(jobs) do (fi, fam, n, inst)
        edges = FAMILIES[fi][2](MersenneTwister(instance_seed(fi, n, inst)), n)
        m = length(edges); costs = maxcut_costs(n, edges); copt = maximum(costs)
        N = get_homogeneous_distribution_from_costs_direct(costs, m, n)
        counts = zeros(Int, m + 1); for c in costs; counts[Int(c)+1] += 1; end
        (; fam, n, inst, m, costs, copt, N, counts)
    end
    println("prep done ($(length(prep)) instances)")
    out = Vector{Vector{String}}(undef, length(prep))
    @threads for k in eachindex(prep)
        q = prep[k]
        out[k] = instance_rows(q.fam, q.n, q.inst, q.m, q.costs, q.copt, q.N, q.counts)
        println("done: ", q.fam, " n=", q.n, " inst=", q.inst)
    end
    path = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,m,p,regret,e_star,e_hat,delta_e,margin,cancel_factor")
        for rs in out, r in rs; println(io, r); end
    end
    println("wrote $path; decomposition identity held everywhere")
end

main()
