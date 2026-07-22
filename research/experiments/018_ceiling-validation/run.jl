#=
E018 — Is the grid ceiling a fair reference for regret?

The journal-readiness review flagged two objections to the regret metric:
(a) the p=1 40x40 grid and especially the p=3 8^4 linear-ramp grid may sit
below the true optimum, understating every regret; (b) no external baseline
(parameter transfer) is compared. This experiment measures both.

For each instance (same families/sizes/seeds as E002/E004/E014):
  1. p=1: grid ceiling, then Nelder-Mead refinement of the true objective
     from the grid argmax (2 params). gap1 = refined - grid.
  2. p=3: 8^4 ramp-grid ceiling, then NM over the 4 ramp endpoints.
     gap3_ramp = refined - grid.
  3. p=3: NM over all 6 angles (gamma_1..3, beta_1..3) from the refined
     ramp, testing the linear-ramp restriction itself.
     gap3_full = refined6 - refined4.
  4. Transfer baseline: the mean true-argmax angles of the ER(0.5), n=12
     source cell, applied verbatim to every instance (Galda-style small-to-
     large/cross-family transfer). ar_transfer_p1/p3.

Run from the repo root (threaded, ~10-20 min):
  JULIA_NUM_THREADS=auto julia --project research/experiments/018_ceiling-validation/run.jl
Smoke test: E18_SMOKE=1 julia --project .../run.jl
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E18_SMOKE", "0") == "1"
const SEED  = 20260611                       # identical instances to exp 002/004/014
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
const P1_sched  = vec([([g], [b]) for g in P1γ, b in P1β])
const R_combos  = vec([(g1,gf,b1,bf) for g1 in Rγ, gf in Rγ, b1 in Rβ, bf in Rβ])
const R_sched   = [linear_ramp(c..., 3) for c in R_combos]

# ---------------- Nelder-Mead (box-clamped, minimizes f) ----------------
function nelder_mead(f, x0::Vector{Float64}, lo::Vector{Float64}, hi::Vector{Float64};
                     step=0.05, maxit=400, ftol=1e-9)
    clamp!(v) = (for i in eachindex(v); v[i] = min(max(v[i], lo[i]), hi[i]); end; v)
    d = length(x0)
    simplex = [copy(x0)]
    for i in 1:d
        v = copy(x0)
        v[i] += step * (hi[i] - lo[i]) * (v[i] + step > hi[i] ? -1 : 1)
        push!(simplex, clamp!(v))
    end
    fv = [f(v) for v in simplex]
    for _ in 1:maxit
        ord = sortperm(fv); simplex = simplex[ord]; fv = fv[ord]
        (fv[end] - fv[1]) < ftol && break
        centroid = sum(simplex[1:d]) / d
        xr = clamp!(centroid .+ (centroid .- simplex[end]))
        fr = f(xr)
        if fr < fv[1]
            xe = clamp!(centroid .+ 2 .* (centroid .- simplex[end]))
            fe = f(xe)
            (simplex[end], fv[end]) = fe < fr ? (xe, fe) : (xr, fr)
        elseif fr < fv[d]
            simplex[end], fv[end] = xr, fr
        else
            xc = clamp!(centroid .+ 0.5 .* (simplex[end] .- centroid))
            fc = f(xc)
            if fc < fv[end]
                simplex[end], fv[end] = xc, fc
            else
                for i in 2:d+1
                    simplex[i] = clamp!(simplex[1] .+ 0.5 .* (simplex[i] .- simplex[1]))
                    fv[i] = f(simplex[i])
                end
            end
        end
    end
    ib = argmin(fv)
    return simplex[ib], fv[ib]
end

# ---------------- per-instance work ----------------
function instance_row(fam, n, inst, m, costs, copt, t1, t3)
    ar(γs, βs) = qaoa_expectation(costs, n, γs, βs) / copt

    # p=1 grid + refinement
    L1 = [ar(s...) for s in P1_sched]
    c1 = maximum(L1); i1 = argmax(L1)
    g0, b0 = P1_sched[i1][1][1], P1_sched[i1][2][1]
    x1, f1 = nelder_mead(x -> -ar([x[1]], [x[2]]), [g0, b0], [0.0, 0.0], [Float64(π), π/2])
    r1 = -f1

    # p=3 ramp grid + 4-param refinement
    L3 = [ar(s...) for s in R_sched]
    c3 = maximum(L3); i3 = argmax(L3)
    c0 = collect(Float64, R_combos[i3])
    x3, f3 = nelder_mead(x -> -ar(linear_ramp(x..., 3)...), c0,
                         [0.05, 0.05, 0.05, 0.05], [1.6, 1.6, 0.8, 0.8])
    r3 = -f3

    # p=3 full 6-parameter refinement from the refined ramp
    γr, βr = linear_ramp(x3..., 3)
    x6, f6 = nelder_mead(x -> -ar(x[1:3], x[4:6]), vcat(γr, βr),
                         fill(0.0, 6), vcat(fill(Float64(π), 3), fill(Float64(π)/2, 3)))
    r6 = -f6

    # transfer baseline
    art1 = ar([t1[1]], [t1[2]])
    art3 = ar(linear_ramp(t3..., 3)...)

    return join((fam, n, inst, m,
        round(c1, digits=6), round(r1, digits=6), round(r1 - c1, digits=6),
        round(c3, digits=6), round(r3, digits=6), round(r3 - c3, digits=6),
        round(r6, digits=6), round(r6 - r3, digits=6),
        round(art1, digits=6), round(art3, digits=6)), ",")
end

function main()
    jobs = [(fi, fam, n, inst) for (fi, (fam, _)) in enumerate(FAMILIES)
            for n in NS for inst in 1:INST]
    prep = map(jobs) do (fi, fam, n, inst)
        edges = FAMILIES[fi][2](MersenneTwister(instance_seed(fi, n, inst)), n)
        costs = maxcut_costs(n, edges)
        (; fam, n, inst, m = length(edges), costs, copt = maximum(costs))
    end

    # transfer source: true argmax angles of the ER(0.5), n=NS[1] cell, averaged
    src = filter(q -> q.fam == "ER(0.5)" && q.n == NS[1], prep)
    isempty(src) && error("no transfer-source instances found")
    t1s = Vector{Vector{Float64}}(undef, length(src))
    t3s = Vector{Vector{Float64}}(undef, length(src))
    @threads for k in eachindex(src)
        q = src[k]
        ar(γs, βs) = qaoa_expectation(q.costs, q.n, γs, βs) / q.copt
        i1 = argmax([ar(s...) for s in P1_sched])
        i3 = argmax([ar(s...) for s in R_sched])
        t1s[k] = [P1_sched[i1][1][1], P1_sched[i1][2][1]]
        t3s[k] = collect(Float64, R_combos[i3])
    end
    t1 = sum(t1s) / length(t1s)
    t3 = sum(t3s) / length(t3s)
    println("transfer angles p1 (γ,β) = ", round.(t1, digits=4),
            "  p3 ramp = ", round.(t3, digits=4))

    rows = Vector{String}(undef, length(prep))
    @threads for k in eachindex(prep)
        q = prep[k]
        rows[k] = instance_row(q.fam, q.n, q.inst, q.m, q.costs, q.copt, t1, t3)
        println("done: ", q.fam, " n=", q.n, " inst=", q.inst)
    end
    path = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,m,ceil_p1,refined_p1,gap_p1," *
                    "ceil_p3,refined_p3_ramp,gap_p3_ramp," *
                    "refined_p3_full,gap_p3_full,ar_transfer_p1,ar_transfer_p3")
        for r in rows; println(io, r); end
    end
    println("wrote $path ($(length(rows)) instances)")
end

main()
