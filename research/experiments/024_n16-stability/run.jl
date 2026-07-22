#=
E024 — Do the headline parameter-space findings persist at n = 16?

E014/E017 (argmax displacement predicts regret), E020 (normalization halves
depth regret), and E022 (signed anatomy: uniform overprediction, winner's
curse, margin absorption) were established at n = 12, 14. E010 extended
REGRET to n = 16-18 but not the landscape geometry. This experiment
re-measures, at n = 16 (7 families x 5 instances, p = 1 and 3):

  - regret (normalized and unnormalized proxy objective),
  - argmax displacement and fidelity deficit at the chosen angles,
  - the signed decomposition e(theta*), e(theta-hat), margin.

Analysis (analyze.jl or the README) then reports: Spearman(regret, argmax
displacement) and (regret, fidelity deficit) at both depths; the
normalization benefit at p = 3; overprediction sign counts; winner's-curse
ratio. Instances use the same seed scheme as E002/E010.

Run: JULIA_NUM_THREADS=auto julia --project research/experiments/024_n16-stability/run.jl  (long: ~1-2 h)
Smoke: E24_SMOKE=1 ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E24_SMOKE", "0") == "1"
const SEED  = 20260611
const NVAL  = SMOKE ? 10 : 16
const INST  = SMOKE ? 1 : 5
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
const P1_coord = vec([(g/π, b/(π/2)) for g in P1γ, b in P1β])
const R_combos = vec([(g1,gf,b1,bf) for g1 in Rγ, gf in Rγ, b1 in Rβ, bf in Rβ])
const R_sched  = [linear_ramp(c..., 3) for c in R_combos]
const R_coord  = [ ((c[1]-0.05)/1.55, (c[2]-0.05)/1.55, (c[3]-0.05)/0.75, (c[4]-0.05)/0.75)
                   for c in R_combos ]
normdist(a::Tuple, b::Tuple) = sqrt(sum((a .- b) .^ 2)) / sqrt(length(a))

function proxy_cols(N, γmat, βmat, counts, copt)
    Q = QAOA_proxy_multi(N, γmat, βmat)[end]
    K = size(Q, 2)
    raw = zeros(Float64, K); nrm = zeros(Float64, K)
    for k in 1:K, ci in axes(Q, 1)
        w = counts[ci] * abs2(Q[ci, k])
        raw[k] += w * (ci - 1); nrm[k] += w
    end
    return raw ./ nrm ./ copt, raw ./ copt   # normalized, unnormalized objectives
end

function instance_rows(fam, n, inst, m, costs, copt, N, counts)
    rows = String[]
    for (p, scheds, coords, γmat, βmat) in (
            (1, P1_sched, P1_coord,
                reshape([s[1][1] for s in P1_sched], :, 1),
                reshape([s[2][1] for s in P1_sched], :, 1)),
            (3, R_sched, R_coord,
                linear_ramp_matrix([c[1] for c in R_combos], [c[2] for c in R_combos],
                                   [c[3] for c in R_combos], [c[4] for c in R_combos], 3)...))
        F = [qaoa_expectation(costs, n, γs, βs) / copt for (γs, βs) in scheds]
        fh_norm, fh_raw = proxy_cols(N, γmat, βmat, counts, copt)
        istar = argmax(F)
        ihat_n = argmax(fh_norm); ihat_r = argmax(fh_raw)
        regret_n = F[istar] - F[ihat_n]
        regret_r = F[istar] - F[ihat_r]
        disp = normdist(coords[ihat_n], coords[istar])
        γe, βe = scheds[ihat_n]
        ov = compressed_qaoa_trajectory(costs, n, γe, βe).overlap[end]
        estar = F[istar] - fh_norm[istar]
        ehat  = F[ihat_n] - fh_norm[ihat_n]
        margin = fh_norm[ihat_n] - fh_norm[istar]
        push!(rows, join((fam, n, inst, m, p,
            round(F[istar], digits=6), round(regret_n, digits=6),
            round(regret_r, digits=6), round(disp, digits=6),
            round(1 - ov, digits=6), round(estar, digits=6),
            round(ehat, digits=6), round(margin, digits=6)), ","))
    end
    return rows
end

function main()
    jobs = [(fi, fam, inst) for (fi, (fam, _)) in enumerate(FAMILIES) for inst in 1:INST]
    prep = map(jobs) do (fi, fam, inst)
        edges = FAMILIES[fi][2](MersenneTwister(instance_seed(fi, NVAL, inst)), NVAL)
        m = length(edges); costs = maxcut_costs(NVAL, edges); copt = maximum(costs)
        N = get_homogeneous_distribution_from_costs_direct(costs, m, NVAL)
        counts = zeros(Int, m + 1); for c in costs; counts[Int(c)+1] += 1; end
        (; fam, inst, m, costs, copt, N, counts)
    end
    println("prep done ($(length(prep)) instances at n=$NVAL)")
    out = Vector{Vector{String}}(undef, length(prep))
    @threads for k in eachindex(prep)
        q = prep[k]
        out[k] = instance_rows(q.fam, NVAL, q.inst, q.m, q.costs, q.copt, q.N, q.counts)
        println("done: ", q.fam, " inst=", q.inst)
    end
    path = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,m,p,ceiling,regret_norm,regret_raw," *
                    "argmax_disp,fidelity_deficit,e_star,e_hat,margin")
        for rs in out, r in rs; println(io, r); end
    end
    println("wrote $path")
end

main()
