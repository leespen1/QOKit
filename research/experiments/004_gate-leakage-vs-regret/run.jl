#=
E1.3 (GATE) — Do leakage, fidelity, and parameter-setting regret rank graph
families identically, and how much regret does the analytical N add on top?

For every instance (same seeds as experiment 002), at p=1 (full grid) and p=3
(linear-ramp endpoint grid):
  - ceiling: best real-QAOA AR over the grid;
  - exact-compression proxy (same-instance empirical N and P): argmax → real AR
    → regret_emp (pure compression error in the chosen-parameter sense);
  - analytical PaperProxy N (Binomial/Multinomial with p_eff = 2m/(n(n-1))) and
    binomial P: argmax → real AR → regret_paper (compression + model error);
  - leakage diagnostics along the empirical proxy's chosen schedule:
    Σλ, final overlap and distance (compressed_qaoa_trajectory).

Go/no-go (Phase-1 gate): families ranked by mean Σλ should rank by mean
regret_emp with |Spearman ρ| ≳ 0.8. The regret_paper − regret_emp gap measures
model error and tests whether non-ER failure is a model-error story.

Run from the repo root (~1 h with threads):
  JULIA_NUM_THREADS=auto julia --project research/experiments/004_gate-leakage-vs-regret/run.jl
Smoke test: E1_SMOKE=1 julia --project ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads
using Statistics: mean

const SMOKE = get(ENV, "E1_SMOKE", "0") == "1"

const SEED = 20260611   # shared with experiment 002 → identical instances
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

instance_seed(fam_idx, n, inst) = SEED + 10_000 * fam_idx + 100 * n + inst

"Best real-QAOA expectation over schedules (serial — outer loop is threaded)."
function grid_ceiling(costs, n, schedules)
    best = -Inf
    for (γs, βs) in schedules
        best = max(best, qaoa_expectation(costs, n, γs, βs))
    end
    return best
end

"Argmax index of a proxy over K schedules given as K×p matrices."
function proxy_argmax(N, P, n, γmat, βmat)
    Qs = QAOA_proxy_multi(N, γmat, βmat)
    return argmax(vec(expectation(Qs[end], P, n)))
end

function main()
    # Enumerate all jobs and pre-generate edges (cheap) so PaperProxy N arrays
    # can be cached per distinct (n, m) — they are expensive and shared widely.
    jobs = []
    for (fam_idx, (fam_name, gen)) in enumerate(FAMILIES), n in NS, inst in 1:INSTANCES
        seed = instance_seed(fam_idx, n, inst)
        edges = gen(MersenneTwister(seed), n)
        push!(jobs, (; fam_name, n, inst, seed, edges))
    end

    nm_pairs = unique((j.n, length(j.edges)) for j in jobs)
    println("distinct (n, m) pairs for PaperProxy N: ", length(nm_pairs))
    paper_cache = Dict{Tuple{Int,Int},Tuple{Array{Float64,3},Vector{Float64}}}()
    cache_lock = ReentrantLock()
    @threads for (n, m) in nm_pairs
        p_eff = 2m / (n * (n - 1))
        proxy = PaperProxy(m, n, p_eff)
        N_paper = cpu_compute_homodist(proxy)
        P_paper = [P_cost_distribution(proxy, c) for c in 0:m]
        lock(() -> (paper_cache[(n, m)] = (N_paper, P_paper)), cache_lock)
    end
    println("PaperProxy cache built")

    # Schedule grids (shared across instances)
    p1_schedules = vec([([γ], [β]) for γ in P1_γ, β in P1_β])
    γmat1 = reshape([s[1][1] for s in p1_schedules], :, 1)
    βmat1 = reshape([s[2][1] for s in p1_schedules], :, 1)
    ramp_combos = vec([(g1, gf, b1, bf) for g1 in RAMP_γ, gf in RAMP_γ,
                                            b1 in RAMP_β, bf in RAMP_β])
    ramp_schedules = [linear_ramp(c..., P_RAMP) for c in ramp_combos]
    γmat3, βmat3 = linear_ramp_matrix(
        [c[1] for c in ramp_combos], [c[2] for c in ramp_combos],
        [c[3] for c in ramp_combos], [c[4] for c in ramp_combos], P_RAMP)

    # Empirical N arrays serially up front: the homodist function is internally
    # threaded (@threads :static), which cannot nest inside the threaded loop below.
    homodists = Vector{Array{Float64,3}}(undef, length(jobs))
    all_costs = Vector{Vector{Float64}}(undef, length(jobs))
    for k in eachindex(jobs)
        j = jobs[k]
        all_costs[k] = maxcut_costs(j.n, j.edges)
        homodists[k] = get_homogeneous_distribution_from_costs_direct(
            all_costs[k], length(j.edges), j.n)
    end
    println("empirical homodists built")

    results = Vector{String}(undef, length(jobs))
    @threads for k in eachindex(jobs)
        j = jobs[k]
        n, m = j.n, length(j.edges)
        costs = all_costs[k]
        c_opt = maximum(costs)

        N_emp = homodists[k]
        counts = zeros(Int, m + 1)
        for cst in costs
            counts[Int(cst) + 1] += 1
        end
        P_emp = counts ./ (1 << n)
        N_paper, P_paper = paper_cache[(n, m)]

        row = Any[j.fam_name, n, j.inst, j.seed, m]
        for (γmat, βmat, schedules, p) in ((γmat1, βmat1, p1_schedules, 1),
                                           (γmat3, βmat3, ramp_schedules, P_RAMP))
            ceiling = grid_ceiling(costs, n, schedules) / c_opt

            i_emp = proxy_argmax(N_emp, P_emp, n, γmat, βmat)
            γ_emp, β_emp = schedules[i_emp]
            ar_emp = qaoa_expectation(costs, n, γ_emp, β_emp) / c_opt

            i_pap = proxy_argmax(N_paper, P_paper, n, γmat, βmat)
            γ_pap, β_pap = schedules[i_pap]
            ar_pap = qaoa_expectation(costs, n, γ_pap, β_pap) / c_opt

            traj = compressed_qaoa_trajectory(costs, n, γ_emp, β_emp)
            row = vcat(row, [ceiling, ar_emp, ar_pap,
                             sum(traj.leakage), traj.overlap[end], traj.distance[end]])
        end
        results[k] = join(row, ",")
        println("done: $(j.fam_name) n=$n inst=$(j.inst)")
    end

    outpath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m," *
                    "ceil_p1,ar_emp_p1,ar_paper_p1,sum_leakage_p1,overlap_p1,distance_p1," *
                    "ceil_p3,ar_emp_p3,ar_paper_p3,sum_leakage_p3,overlap_p3,distance_p3")
        foreach(r -> println(io, r), results)
    end
    println("\nE1.3 sweep complete → $outpath")
end

main()
