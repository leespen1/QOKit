#=
E1.3 follow-up #2 — Does rejecting only *impossible* predictions repair the
analytical proxy on dense ER(0.5) without harming the sparse families?

Experiment 005 showed norm-tolerance filters are catastrophic off ER(0.5):
the analytical model's calibration (values, norms) is broken on sparse graphs
even where its argmax is nearly perfect. The surviving idea: a prediction
⟨C⟩ > m exceeds the largest possible cut and is *certainly* an artifact, so
restrict the argmax to grid points with predicted ⟨C⟩ ≤ m (and, as a second
variant, ≤ m/2 + headroom·m/2 — between the random-guess value and the
physical maximum). Unlike norm tolerance, this never rejects a well-located
peak with miscalibrated height unless the height is impossible.

Run from the repo root (~30 min with threads):
  JULIA_NUM_THREADS=auto julia --project research/experiments/006_physicality-filter/run.jl
Smoke test: E1_SMOKE=1 julia --project ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E1_SMOKE", "0") == "1"

const SEED = 20260611   # shared with experiments 002/004/005 → identical instances
const NS = SMOKE ? [10] : [12, 14]
const INSTANCES = SMOKE ? 2 : 30
const P1_GRID_LEN = SMOKE ? 12 : 40
const RAMP_GRID_LEN = SMOKE ? 4 : 8
const P_RAMP = 3

const P1_γ = collect(range(0.0, π; length=P1_GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=P1_GRID_LEN))
const RAMP_γ = collect(range(0.05, 1.6; length=RAMP_GRID_LEN))
const RAMP_β = collect(range(0.05, 0.8; length=RAMP_GRID_LEN))

const FAMILIES = [
    ("ER(0.5)",       (rng, n) -> erdos_renyi_edges(n, 0.5; rng)),
    ("ER(0.25)",      (rng, n) -> erdos_renyi_edges(n, 0.25; rng)),
    ("BA(k=2)",       (rng, n) -> barabasi_albert_edges(n, 2; rng)),
    ("BA(k=4)",       (rng, n) -> barabasi_albert_edges(n, 4; rng)),
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("WS(k=4;b=0.5)", (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]

function grid_ceiling(costs, n, schedules)
    best = -Inf
    for (γs, βs) in schedules
        best = max(best, qaoa_expectation(costs, n, γs, βs))
    end
    return best
end

function main()
    jobs = []
    for (fam_idx, (fam_name, gen)) in enumerate(FAMILIES), n in NS, inst in 1:INSTANCES
        seed = SEED + 10_000 * fam_idx + 100 * n + inst
        edges = gen(MersenneTwister(seed), n)
        push!(jobs, (; fam_name, n, inst, seed, edges))
    end

    nm_pairs = unique((j.n, length(j.edges)) for j in jobs)
    paper_cache = Dict{Tuple{Int,Int},Tuple{Array{Float64,3},Vector{Float64}}}()
    cache_lock = ReentrantLock()
    @threads for (n, m) in nm_pairs
        proxy = PaperProxy(m, n, 2m / (n * (n - 1)))
        N_paper = cpu_compute_homodist(proxy)
        P_paper = [P_cost_distribution(proxy, c) for c in 0:m]
        lock(() -> (paper_cache[(n, m)] = (N_paper, P_paper)), cache_lock)
    end
    println("PaperProxy cache built (", length(nm_pairs), " (n, m) pairs)")

    p1_schedules = vec([([γ], [β]) for γ in P1_γ, β in P1_β])
    γmat1 = reshape([s[1][1] for s in p1_schedules], :, 1)
    βmat1 = reshape([s[2][1] for s in p1_schedules], :, 1)
    ramp_combos = vec([(g1, gf, b1, bf) for g1 in RAMP_γ, gf in RAMP_γ,
                                            b1 in RAMP_β, bf in RAMP_β])
    ramp_schedules = [linear_ramp(c..., P_RAMP) for c in ramp_combos]
    γmat3, βmat3 = linear_ramp_matrix(
        [c[1] for c in ramp_combos], [c[2] for c in ramp_combos],
        [c[3] for c in ramp_combos], [c[4] for c in ramp_combos], P_RAMP)

    results = Vector{String}(undef, length(jobs))
    @threads for k in eachindex(jobs)
        j = jobs[k]
        n, m = j.n, length(j.edges)
        costs = maxcut_costs(n, j.edges)
        c_opt = maximum(costs)
        N_paper, P_paper = paper_cache[(n, m)]

        row = Any[j.fam_name, n, j.inst, j.seed, m]
        for (γmat, βmat, schedules) in ((γmat1, βmat1, p1_schedules),
                                        (γmat3, βmat3, ramp_schedules))
            ceiling = grid_ceiling(costs, n, schedules) / c_opt
            vals = vec(expectation(QAOA_proxy_multi(N_paper, γmat, βmat)[end], P_paper, n))

            # When a cap rejects every grid point the model gives no usable
            # signal; fall back to the raw argmax and let frac_* record it.
            ar_at = function (mask)
                masked = copy(vals)
                any(mask) && (masked[.!mask] .= -Inf)
                return qaoa_expectation(costs, n, schedules[argmax(masked)]...) / c_opt
            end
            ar_raw = ar_at(trues(length(vals)))
            # Physical cap: no cut exceeds m edges, so vals > m are certain artifacts
            ar_phys = ar_at(vals .<= m)
            # Sensitivity check with a tighter (heuristic) cap; can over-reject on
            # sparse graphs whose optima approach m
            ar_tight = ar_at(vals .<= 0.75 * m)
            frac_phys = count(vals .<= m) / length(vals)
            frac_tight = count(vals .<= 0.75 * m) / length(vals)

            row = vcat(row, [ceiling, ar_raw, ar_phys, ar_tight, frac_phys, frac_tight])
        end
        results[k] = join(row, ",")
        println("done: $(j.fam_name) n=$n inst=$(j.inst)")
    end

    outpath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m," *
                    "ceil_p1,ar_raw_p1,ar_phys_p1,ar_tight_p1,frac_phys_p1,frac_tight_p1," *
                    "ceil_p3,ar_raw_p3,ar_phys_p3,ar_tight_p3,frac_phys_p3,frac_tight_p3")
        foreach(r -> println(io, r), results)
    end
    println("\nE1.3-followup-2 sweep complete → $outpath")
end

main()
