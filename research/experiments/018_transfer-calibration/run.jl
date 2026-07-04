#=
E018 — Transfer calibration (paper-critique T2.4): how does proxy regret
compare to the standard cheap alternative, transferring angles optimized on a
SMALL instance of the same family?

For every target instance already measured in exps 002 (n=12,14) and 010
(n=16,18) — regenerated exactly from the seeds recorded in their CSVs, with
the edge count asserted against the CSV — evaluate real QAOA at angles
grid-optimized on the TRUE landscape of a small n=10 source instance of the
same family (3 independent sources per family), at p=1 and p=3 ramps, on the
same schedule grids as 002/010. Regret is against the target's own recorded
grid ceiling, so transfer and proxy numbers are directly comparable per cell.

Output: results.csv, long format — one row per (target instance, source, p).

Submit:  cd research/experiments/018_transfer-calibration && sbatch run.sb
Smoke:   E18_SMOKE=1 julia --project research/experiments/018_transfer-calibration/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E18_SMOKE", "0") == "1"

const SOURCE_SEED = 20260704
const SOURCE_N = 10
const N_SOURCES = SMOKE ? 1 : 3
const GRID_LEN = 40
const RAMP_LEN = 8
const P_RAMP = 3

const P1_γ = collect(range(0.0, π; length=GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=GRID_LEN))
const RAMP_γ = collect(range(0.05, 1.6; length=RAMP_LEN))
const RAMP_β = collect(range(0.05, 0.8; length=RAMP_LEN))

const FAMILIES = Dict(
    "ER(0.5)"       => (rng, n) -> erdos_renyi_edges(n, 0.5; rng),
    "ER(0.25)"      => (rng, n) -> erdos_renyi_edges(n, 0.25; rng),
    "BA(k=2)"       => (rng, n) -> barabasi_albert_edges(n, 2; rng),
    "BA(k=4)"       => (rng, n) -> barabasi_albert_edges(n, 4; rng),
    "WS(k=4;b=0.1)" => (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng),
    "WS(k=4;b=0.5)" => (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng),
    "3-regular"     => (rng, n) -> random_regular_edges(n, 3; rng),
)

struct Target
    family::String
    n::Int
    instance::Int
    seed::Int
    m::Int
    ceil_p1::Float64
    ceil_p3::Float64
end

"Read targets from exp 002 (cols 5,9,13 = m, ceil_p1, ceil_p3) and exp 010 (cols 5,8,11)."
function load_targets()
    targets = Target[]
    f002 = joinpath(@__DIR__, "..", "002_baselines-and-headroom", "results.csv")
    for (i, line) in enumerate(eachline(f002))
        i == 1 && continue
        p = split(line, ",")
        push!(targets, Target(p[1], parse(Int, p[2]), parse(Int, p[3]),
                              parse(Int, p[4]), parse(Int, p[5]),
                              parse(Float64, p[9]), parse(Float64, p[13])))
    end
    dir010 = joinpath(@__DIR__, "..", "010_scaleup-ranking")
    for f in sort(filter(f -> occursin(r"^results_task\d+\.csv$", f), readdir(dir010)))
        for (i, line) in enumerate(eachline(joinpath(dir010, f)))
            i == 1 && continue
            p = split(line, ",")
            push!(targets, Target(p[1], parse(Int, p[2]), parse(Int, p[3]),
                                  parse(Int, p[4]), parse(Int, p[5]),
                                  parse(Float64, p[8]), parse(Float64, p[11])))
        end
    end
    targets
end

"Grid-optimal (γ,β) at p=1 and ramp endpoints at p=3 on one instance's TRUE landscape."
function source_optimal_angles(costs, n)
    best1 = (-Inf, 0.0, 0.0)
    for γ in P1_γ, β in P1_β
        v = qaoa_expectation(costs, n, [γ], [β])
        v > best1[1] && (best1 = (v, γ, β))
    end
    best3 = (-Inf, 0.0, 0.0, 0.0, 0.0)
    for g1 in RAMP_γ, gf in RAMP_γ, b1 in RAMP_β, bf in RAMP_β
        γs, βs = linear_ramp(g1, gf, b1, bf, P_RAMP)
        v = qaoa_expectation(costs, n, γs, βs)
        v > best3[1] && (best3 = (v, g1, gf, b1, bf))
    end
    return (best1[2], best1[3]), (best3[2], best3[3], best3[4], best3[5])
end

function main()
    targets = load_targets()
    SMOKE && (targets = filter(t -> t.instance <= 2 && t.n <= 12, targets))
    println("targets: ", length(targets))

    # per-family sources, optimized once
    sources = Dict{String, Vector{Tuple{Tuple, Tuple}}}()
    for fam in keys(FAMILIES)
        reps = Vector{Tuple{Tuple, Tuple}}()
        for rep in 1:N_SOURCES
            seed = SOURCE_SEED + 1000 * findfirst(==(fam), sort(collect(keys(FAMILIES)))) + rep
            edges = FAMILIES[fam](MersenneTwister(seed), SOURCE_N)
            costs = maxcut_costs(SOURCE_N, edges)
            push!(reps, source_optimal_angles(costs, SOURCE_N))
        end
        sources[fam] = reps
        println("sources ready: $fam")
        flush(stdout)
    end

    rows = Vector{String}(undef, length(targets) * N_SOURCES * 2)
    @threads for k in eachindex(targets)
        t = targets[k]
        edges = FAMILIES[t.family](MersenneTwister(t.seed), t.n)
        @assert length(edges) == t.m "seed regeneration mismatch: $(t.family) n=$(t.n) inst=$(t.instance) m=$(length(edges)) vs CSV $(t.m)"
        costs = maxcut_costs(t.n, edges)
        c_opt = maximum(costs)
        for (rep, ((γ1, β1), (g1, gf, b1, bf))) in enumerate(sources[t.family])
            ar1 = qaoa_expectation(costs, t.n, [γ1], [β1]) / c_opt
            γs, βs = linear_ramp(g1, gf, b1, bf, P_RAMP)
            ar3 = qaoa_expectation(costs, t.n, γs, βs) / c_opt
            @assert t.ceil_p1 - ar1 > -1e-8 && t.ceil_p3 - ar3 > -1e-8 "transfer beat the recorded ceiling"
            base = (k - 1) * N_SOURCES * 2 + (rep - 1) * 2
            rows[base + 1] = join(Any[t.family, t.n, t.instance, t.seed, rep, 1,
                                      ar1, t.ceil_p1, t.ceil_p1 - ar1], ",")
            rows[base + 2] = join(Any[t.family, t.n, t.instance, t.seed, rep, 3,
                                      ar3, t.ceil_p3, t.ceil_p3 - ar3], ",")
        end
        k % 50 == 0 && (println("targets done: $k"); flush(stdout))
    end

    outpath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,source_rep,p,ar_transfer,ar_ceiling,regret_transfer")
        foreach(r -> println(io, r), rows)
    end
    println("E018 complete → $outpath")
end

main()
