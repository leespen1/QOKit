#=
E039 — Is the proxy's overprediction predictable, and does correcting it
help parameter setting?

E038 found the normalized proxy objective OVERPREDICTS the true objective
essentially everywhere, most at its own argmax (winner's curse). If the
signed error e(theta) = F - F-hat is predictable from a proxy-side
observable, subtracting the prediction gives a corrected objective at zero
extra cost. The natural candidate is the norm loss 1 - ||phi(theta)||^2
(free from Theorem 2's bookkeeping): leakage is what the compression throws
away, and the normalization redistributes the surviving amplitude.

Design (honest split):
  Stage A (probe): instances inst=1 of each (family, n) cell (14 total),
  p=1 grid. Regress e(theta) on x(theta) = 1 - ||phi||^2 pooled across the
  14 full grids: e ~ a + b x. Report per-instance Pearson correlations and
  the pooled slope b-hat. Only the slope matters for the argmax.
  Stage B (evaluate): corrected objective F-c = F-hat + b-hat * x on ALL
  instances and both depths; report regret(corrected) vs regret(normalized)
  and vs regret(unnormalized), with the probe instances excluded from the
  headline numbers.

Run: JULIA_NUM_THREADS=auto julia --project research/experiments/039_bias-corrected-objective/run.jl  (~12 min)
Smoke: E23_SMOKE=1 ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E23_SMOKE", "0") == "1"
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

"Per-column (normalized objective, norm^2) of a proxy-multi run."
function proxy_cols(N, γmat, βmat, counts, copt)
    Q = QAOA_proxy_multi(N, γmat, βmat)[end]
    K = size(Q, 2)
    fh = zeros(Float64, K); n2 = zeros(Float64, K)
    for k in 1:K
        raw = 0.0; nrm = 0.0
        for ci in axes(Q, 1)
            w = counts[ci] * abs2(Q[ci, k])
            raw += w * (ci - 1); nrm += w
        end
        fh[k] = raw / nrm / copt; n2[k] = nrm
    end
    return fh, n2
end

function prep_instance(fi, fam, n, inst)
    edges = FAMILIES[fi][2](MersenneTwister(instance_seed(fi, n, inst)), n)
    m = length(edges); costs = maxcut_costs(n, edges); copt = maximum(costs)
    N = get_homogeneous_distribution_from_costs_direct(costs, m, n)
    counts = zeros(Int, m + 1); for c in costs; counts[Int(c)+1] += 1; end
    (; fi, fam, n, inst, m, costs, copt, N, counts)
end

pearson(x, y) = begin
    mx = sum(x)/length(x); my = sum(y)/length(y)
    sx = sqrt(sum(abs2, x .- mx)); sy = sqrt(sum(abs2, y .- my))
    sum((x .- mx) .* (y .- my)) / (sx * sy)
end

function main()
    jobs = [(fi, fam, n, inst) for (fi, (fam, _)) in enumerate(FAMILIES)
            for n in NS for inst in 1:INST]
    prep = map(j -> prep_instance(j...), jobs)
    println("prep done ($(length(prep)) instances)")

    # ---------- Stage A: probe regressions on inst == 1, p = 1 ----------
    probes = filter(q -> q.inst == 1, prep)
    γ1 = reshape([s[1][1] for s in P1_sched], :, 1)
    β1 = reshape([s[2][1] for s in P1_sched], :, 1)
    probe_corrs = Vector{Float64}(undef, length(probes))
    XY = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, length(probes))
    @threads for k in eachindex(probes)
        q = probes[k]
        F = [qaoa_expectation(q.costs, q.n, γs, βs) / q.copt for (γs, βs) in P1_sched]
        fh, n2 = proxy_cols(q.N, γ1, β1, q.counts, q.copt)
        e = F .- fh
        x = 1.0 .- n2
        probe_corrs[k] = pearson(e, x)
        XY[k] = (x, e)
    end
    # pooled least-squares slope (with intercept) over all probe grid points
    xs = vcat([xy[1] for xy in XY]...); es = vcat([xy[2] for xy in XY]...)
    mx = sum(xs)/length(xs); me = sum(es)/length(es)
    bhat = sum((xs .- mx) .* (es .- me)) / sum(abs2, xs .- mx)
    println("probe per-instance corr(e, 1-norm^2): ",
            join(string.(round.(probe_corrs, digits=2)), " "))
    println("pooled slope b-hat = ", round(bhat, digits=4))

    # ---------- Stage B: corrected argmax on every instance, both depths ----
    rows = Vector{Vector{String}}(undef, length(prep))
    @threads for k in eachindex(prep)
        q = prep[k]
        out = String[]
        for (p, scheds, γmat, βmat) in (
                (1, P1_sched, γ1, β1),
                (3, R_sched,
                    linear_ramp_matrix([c[1] for c in R_combos], [c[2] for c in R_combos],
                                       [c[3] for c in R_combos], [c[4] for c in R_combos], 3)...))
            F = [qaoa_expectation(q.costs, q.n, γs, βs) / q.copt for (γs, βs) in scheds]
            fh, n2 = proxy_cols(q.N, γmat, βmat, q.counts, q.copt)
            ceilv = maximum(F)
            r_norm = ceilv - F[argmax(fh)]
            fc = fh .+ bhat .* (1.0 .- n2)
            r_corr = ceilv - F[argmax(fc)]
            push!(out, join((q.fam, q.n, q.inst, p, Int(q.inst == 1),
                round(ceilv, digits=6), round(r_norm, digits=6),
                round(r_corr, digits=6)), ","))
        end
        rows[k] = out
        println("done: ", q.fam, " n=", q.n, " inst=", q.inst)
    end
    path = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,p,is_probe,ceiling,regret_normalized,regret_corrected")
        for rs in rows, r in rs; println(io, r); end
    end
    open(joinpath(@__DIR__, SMOKE ? "slope_smoke.txt" : "slope.txt"), "w") do io
        println(io, "bhat,", bhat)
        println(io, "probe_corrs,", join(round.(probe_corrs, digits=4), ","))
    end
    println("wrote $path")
end

main()
