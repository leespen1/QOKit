#=
E1.2 — Does accumulated per-layer leakage tightly track the true-vs-proxy
overlap deficit (Theorem-2 bound tightness)?

Theorem 2 (telescoping): ‖ψ_p − φ_p‖ ≤ Σ_ℓ λ_ℓ, where λ_ℓ is the norm that
leaks out of the cost-class subspace at layer ℓ and φ is the (unnormalized)
exactly-compressed trajectory. Both sides are measured per layer here, for
p = 20 linear-ramp schedules across graph families. If the bound is loose by
>10× systematically, it is true but vacuous and we lean on measured λ directly.

Run from the repo root (~minutes):
  JULIA_NUM_THREADS=auto julia --project research/experiments/003_leakage-vs-overlap/run.jl
Smoke test: E1_SMOKE=1 julia --project ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E1_SMOKE", "0") == "1"

const SEED = 20260612
const NS = SMOKE ? [10] : [12, 14, 16]
const INSTANCES = SMOKE ? 2 : 10
const P = 20

# Linear-ramp schedules (γ₁→γf, β₁→βf), small to aggressively large angles
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
    jobs = [(fam_idx, fam_name, gen, n, inst, ramp)
            for (fam_idx, (fam_name, gen)) in enumerate(FAMILIES)
            for n in NS, inst in 1:INSTANCES, ramp in RAMPS]
    results = Vector{String}(undef, length(jobs))

    @threads for k in eachindex(jobs)
        fam_idx, fam_name, gen, n, inst, ramp = jobs[k]
        # Same instance seeds as experiment 002 where (family, n, inst) overlap
        seed = SEED + 10_000 * fam_idx + 100 * n + inst
        rng = MersenneTwister(seed)
        edges = gen(rng, n)
        costs = maxcut_costs(n, edges)

        γs, βs = linear_ramp(ramp.g1, ramp.gf, ramp.b1, ramp.bf, P)
        traj = compressed_qaoa_trajectory(costs, n, γs, βs)

        cum = cumsum(traj.leakage)
        rows = String[]
        for ℓ in 1:P
            push!(rows, join([fam_name, n, inst, seed, ramp.id,
                              ramp.g1, ramp.gf, ramp.b1, ramp.bf, ℓ,
                              traj.leakage[ℓ], cum[ℓ],
                              traj.distance[ℓ + 1], traj.overlap[ℓ + 1],
                              traj.compressed_norm[ℓ + 1]], ","))
        end
        results[k] = join(rows, "\n")
        println("done: $fam_name n=$n inst=$inst ramp=$(ramp.id)  " *
                "Σλ=$(round(cum[end], digits=4))  ‖ψ−φ‖=$(round(traj.distance[end], digits=4))  " *
                "|⟨ψ|φ⟩|=$(round(traj.overlap[end], digits=4))")
    end

    outpath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,ramp,g1,gf,b1,bf,layer," *
                    "leakage,cum_leakage,distance,overlap,compressed_norm")
        foreach(r -> println(io, r), results)
    end
    println("\nE1.2 sweep complete → $outpath")
end

main()
