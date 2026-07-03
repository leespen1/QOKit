#=
E014 — At p=3 the proxy-chosen schedules equalize leakage across families
(exp 004: leakage⇔regret ρ≈0), so what governs regret there? This experiment
tests the open framing question: is depth regret a *fidelity* phenomenon (how
close the compressed state stays to QAOA) or an *argmax/landscape-robustness*
phenomenon (how far the proxy's argmax sits from the true optimum, and how much
that costs given the flatness of the true landscape)?

It also explains the exp-004 surprise that the analytical (wrong) model can beat
the exact compression off-ER: regret is set by argmax transfer, and a smoother
landscape can transfer its argmax better even with worse fidelity.

For each instance, at p=1 (40×40 grid) and p=3 (8⁴ ramp grid) we recompute the
quantities exp 004 did not save — the TRUE-landscape argmax location and shape:
  ceiling, ar_emp (exact-compression proxy argmax → real AR), regret;
  argmax_disp  = normalized angle distance between proxy argmax and true argmax;
  robust_frac  = fraction of grid within 0.01 AR of the ceiling (peak flatness);
  overlap      = proxy-state fidelity at the chosen angles (the fidelity axis).

Then `analyze.jl` correlates regret against fidelity vs robustness vs argmax
displacement, at p=1 and p=3 separately.

Run from the repo root (threaded; ~5–10 min):
  JULIA_NUM_THREADS=auto julia --project research/experiments/014_argmax-robustness/run.jl
Smoke test: E14_SMOKE=1 julia --project .../run.jl
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads
using Statistics: mean

const SMOKE = get(ENV, "E14_SMOKE", "0") == "1"
const SEED  = 20260611                       # identical instances to exp 002/004
const NS    = SMOKE ? [10] : [12, 14]
const INST  = SMOKE ? 2 : 10
const GLEN  = SMOKE ? 10 : 40
const RLEN  = SMOKE ? 4 : 8
const EPS   = 0.01                           # AR tolerance defining the "flat peak"

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

# p=1 schedules and their (γ,β) endpoint coordinates (normalized to [0,1])
const P1_sched = vec([([g], [b]) for g in P1γ, b in P1β])
# same 2-D comprehension + vec ⇒ identical (column-major) ordering as P1_sched
const P1_coord = vec([(g/π, b/(π/2)) for g in P1γ, b in P1β])
# p=3 ramp schedules and normalized 4-vectors of endpoints
const R_combos = vec([(g1,gf,b1,bf) for g1 in Rγ, gf in Rγ, b1 in Rβ, bf in Rβ])
const R_sched  = [linear_ramp(c..., 3) for c in R_combos]
const R_coord  = [ ((c[1]-0.05)/1.55, (c[2]-0.05)/1.55, (c[3]-0.05)/0.75, (c[4]-0.05)/0.75)
                   for c in R_combos ]

normdist(a::Tuple, b::Tuple) = sqrt(sum((a .- b) .^ 2)) / sqrt(length(a))

function geometry_rows(fam, n, inst, m, costs, copt, N, P)
    rows = String[]
    for (p, scheds, coords, γmat, βmat) in (
            (1, P1_sched, P1_coord,
                reshape([s[1][1] for s in P1_sched], :, 1),
                reshape([s[2][1] for s in P1_sched], :, 1)),
            (3, R_sched, R_coord,
                linear_ramp_matrix([c[1] for c in R_combos], [c[2] for c in R_combos],
                                   [c[3] for c in R_combos], [c[4] for c in R_combos], 3)...))
        Ltrue = [qaoa_expectation(costs, n, γs, βs) / copt for (γs, βs) in scheds]
        ceil = maximum(Ltrue); it = argmax(Ltrue)
        Qs = QAOA_proxy_multi(N, γmat, βmat)
        vproxy = vec(expectation(Qs[end], P, n))
        ie = argmax(vproxy)
        ar_emp = Ltrue[ie]
        regret = ceil - ar_emp
        disp = normdist(coords[ie], coords[it])
        robust = count(>=(ceil - EPS), Ltrue) / length(Ltrue)
        γe, βe = scheds[ie]
        ov = compressed_qaoa_trajectory(costs, n, γe, βe).overlap[end]
        push!(rows, join((fam, n, inst, m, p, round(ceil,digits=6), round(ar_emp,digits=6),
            round(regret,digits=6), round(disp,digits=6), round(robust,digits=6),
            round(ov,digits=6)), ","))
    end
    return rows
end

function main()
    jobs = [(fi, fam, n, inst) for (fi, (fam, _)) in enumerate(FAMILIES)
            for n in NS for inst in 1:INST]
    # serial prep: the homodist builder is internally @threads :static and cannot
    # be nested inside the threaded landscape loop below.
    prep = map(jobs) do (fi, fam, n, inst)
        edges = FAMILIES[fi][2](MersenneTwister(instance_seed(fi, n, inst)), n)
        m = length(edges); costs = maxcut_costs(n, edges); copt = maximum(costs)
        N = get_homogeneous_distribution_from_costs_direct(costs, m, n)
        counts = zeros(Int, m + 1); for c in costs; counts[Int(c)+1] += 1; end
        (; fam, n, inst, m, costs, copt, N, P = counts ./ (1 << n))
    end
    println("prep done ($(length(prep)) instances); computing landscapes…")
    out = Vector{Vector{String}}(undef, length(prep))
    @threads for k in eachindex(prep)
        q = prep[k]
        out[k] = geometry_rows(q.fam, q.n, q.inst, q.m, q.costs, q.copt, q.N, q.P)
        println("done: ", q.fam, " n=", q.n, " inst=", q.inst)
    end
    path = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,m,p,ceiling,ar_emp,regret,argmax_disp,robust_frac,overlap")
        for rs in out, r in rs; println(io, r); end
    end
    println("wrote $path  ($(length(prep)) instances × 2 depths)")
end

main()
