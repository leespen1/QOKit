#=
E3.2 — How does accumulated leakage grow with depth at scale?

For p = 30 linear-ramp schedules at n ∈ {16, 18, 20}: per-layer leakage λ_ℓ,
cumulative Σλ, distance, overlap, and compressed norm of the exact compressed
trajectory (extends experiment 003 to larger n and deeper p; the question is
whether Σλ grows sublinearly in p at working ramps, and how the growth scales
with n and family). Pure CPU — no GPU needed.

Slurm array: tasks 1..3 map to n ∈ {16, 18, 20} (unset/0 = all, for smoke).

Submit:  sbatch research/experiments/011_depth-scaling/run.sb
Smoke:   E1_SMOKE=1 julia --project research/experiments/011_depth-scaling/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E1_SMOKE", "0") == "1"

const SEED = 20260611
const NS = SMOKE ? [12] : [16, 18, 20]
const INSTANCES = SMOKE ? 2 : 10
const P = SMOKE ? 6 : 30

const RAMPS = [
    (id="small",    g1=0.1, gf=0.4, b1=0.4, bf=0.1),
    (id="moderate", g1=0.1, gf=0.8, b1=0.6, bf=0.1),
    (id="large",    g1=0.2, gf=1.6, b1=0.8, bf=0.1),
    (id="extreme",  g1=0.5, gf=2.5, b1=1.2, bf=0.1),
]

const FAMILIES = [
    ("ER(0.5)",       (rng, n) -> erdos_renyi_edges(n, 0.5; rng)),
    ("ER(0.25)",      (rng, n) -> erdos_renyi_edges(n, 0.25; rng)),
    ("BA(k=2)",       (rng, n) -> barabasi_albert_edges(n, 2; rng)),
    ("BA(k=4)",       (rng, n) -> barabasi_albert_edges(n, 4; rng)),
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("WS(k=4;b=0.5)", (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]

function main()
    tid = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
    ns = tid == 0 ? NS : [NS[tid]]
    println("running n = ", ns, " with $(Threads.nthreads()) threads")

    jobs = [(fam_idx, n, inst, ramp)
            for n in ns
            for fam_idx in eachindex(FAMILIES), inst in 1:INSTANCES, ramp in RAMPS]
    results = Vector{String}(undef, length(jobs))

    @threads for k in eachindex(jobs)
        fam_idx, n, inst, ramp = jobs[k]
        fam_name, gen = FAMILIES[fam_idx]
        seed = SEED + 10_000 * fam_idx + 100 * n + inst
        edges = gen(MersenneTwister(seed), n)
        m = length(edges)
        costs = maxcut_costs(n, edges)

        γs, βs = linear_ramp(ramp.g1, ramp.gf, ramp.b1, ramp.bf, P)
        traj = compressed_qaoa_trajectory(costs, n, γs, βs; num_costs=m + 1)

        cum = cumsum(traj.leakage)
        rows = String[]
        for ℓ in 1:P
            push!(rows, join([fam_name, n, inst, seed, m, ramp.id, ℓ,
                              traj.leakage[ℓ], cum[ℓ],
                              traj.distance[ℓ + 1], traj.overlap[ℓ + 1],
                              traj.compressed_norm[ℓ + 1]], ","))
        end
        results[k] = join(rows, "\n")
        println("done: $fam_name n=$n inst=$inst ramp=$(ramp.id)")
        flush(stdout)
    end

    suffix = tid == 0 ? (SMOKE ? "_smoke" : "_all") : "_task$tid"
    outpath = joinpath(@__DIR__, "results$suffix.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,ramp,layer," *
                    "leakage,cum_leakage,distance,overlap,compressed_norm")
        foreach(r -> println(io, r), results)
    end
    println("E3.2 task(s) complete → $outpath")
end

main()
