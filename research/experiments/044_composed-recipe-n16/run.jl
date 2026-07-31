#=
E044 — Does the composed recipe (sampled N + normalized) survive n=16?

E041 established sampled(S=10)+normalized as the best proxy combination at
n in {12,14}; E042/E043 identified the mechanism (in-set sampling noise
disrupting the winner's curse). The paper now carries the claim (Sec. 6.1).
This experiment replicates the four-cell comparison one size up (n=16),
mirroring E040's stability role for the earlier findings, and adds the
S=3 cell that E042 found strongest.

Design: 7 families x n=16 x 5 instances (fixed seeds, E041 scheme);
variants {exact, samp10, samp3} x {raw, norm}; p=1 grid (40x40) and
p=3 ramp grid (8^4); regret vs the true statevector grid ceiling.

Run: JULIA_NUM_THREADS=auto julia --project research/experiments/044_composed-recipe-n16/run.jl (~40 min)
Smoke: E44_SMOKE=1 ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E44_SMOKE", "0") == "1"
const SEED  = 20260611
const N_V   = SMOKE ? 10 : 16
const INST  = SMOKE ? 2 : 5
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
        s10 = sampled_homogeneous_distribution(costs, m, N_V;
                  samples_per_class=10, rng=MersenneTwister(seed + 777))
        s3  = sampled_homogeneous_distribution(costs, m, N_V;
                  samples_per_class=3,  rng=MersenneTwister(seed + 778))
        (; fam, inst, m, costs, copt, counts,
           variants=[("exact", Nex), ("samp10", s10), ("samp3", s3)])
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
            F = [qaoa_expectation(q.costs, N_V, γs, βs) / q.copt for (γs, βs) in scheds]
            ceilv = maximum(F)
            for (kind, N) in q.variants
                fraw, fnorm = objectives(N, γmat, βmat, q.counts, q.copt)
                push!(out, join((q.fam, N_V, q.inst, p, round(ceilv, digits=6), kind,
                    round(ceilv - F[argmax(fraw)], digits=6),
                    round(ceilv - F[argmax(fnorm)], digits=6)), ","))
            end
        end
        rows[k] = out
        println("done: ", q.fam, " inst=", q.inst)
    end
    path = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,p,ceiling,kind,regret_raw,regret_norm")
        for rs in rows, r in rs; println(io, r); end
    end
    println("wrote $path")
end

main()
