#=
E2.3 — Is the QAOA trajectory low-rank, and is the cost-class subspace the
right low-dimensional subspace for it?

The trajectory ψ_0 … ψ_p spans at most p+1 dimensions, so it is *always*
"low-rank" relative to 2^n. The discriminating questions are:
  1. What fraction E_cc of the trajectory's energy lies inside the
     (m+1)-dimensional cost-class subspace?
  2. How many trajectory-PCA dimensions k_match achieve that same captured
     energy? If k_match ≪ m+1 while E_cc < 1, the cost-class subspace is the
     *wrong* subspace (mis-aimed), not too small.
  3. What is the trajectory's effective dimension (k90/k99 = PCA dims for
     90%/99% energy)?

This separates two stories for where the proxy fails at depth/large angles:
"the state genuinely needs many dimensions" vs "few dimensions suffice but
cost classes aren't them."

Run from the repo root (~10 min with threads):
  JULIA_NUM_THREADS=auto julia --project research/experiments/009_trajectory-pca/run.jl
Smoke test: E1_SMOKE=1 julia --project ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads
using LinearAlgebra: svdvals, svd, norm

const SMOKE = get(ENV, "E1_SMOKE", "0") == "1"

const SEED = 20260611
const NS = SMOKE ? [10] : [12, 14]
const INSTANCES = SMOKE ? 2 : 10
const P = 20

const RAMPS = [
    (id="small",    g1=0.1, gf=0.4, b1=0.4, bf=0.1),
    (id="moderate", g1=0.1, gf=0.8, b1=0.6, bf=0.1),
    (id="large",    g1=0.2, gf=1.6, b1=0.8, bf=0.1),
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
        seed = SEED + 10_000 * fam_idx + 100 * n + inst
        edges = gen(MersenneTwister(seed), n)
        m = length(edges)
        costs = maxcut_costs(n, edges)

        γs, βs = linear_ramp(ramp.g1, ramp.gf, ramp.b1, ramp.bf, P)
        states = qaoa_statevector(costs, n, γs, βs; return_intermediates=true)
        A = stack(states)                       # 2^n × (p+1), unit columns

        # Energy captured by the cost-class subspace
        total = Float64(size(A, 2))             # Σ‖ψ_ℓ‖² = p+1
        captured_cc = 0.0
        for ψ in states
            proj = project_onto_cost_classes(ψ, costs; num_costs=m + 1)
            captured_cc += 1.0 - proj.residual_norm^2
        end
        E_cc = captured_cc / total

        # Trajectory PCA: optimal rank-k captured energy from singular values
        σ2 = svdvals(A) .^ 2
        cum = cumsum(σ2) ./ sum(σ2)
        k_match = something(findfirst(>=(E_cc), cum), length(cum))
        k90 = something(findfirst(>=(0.90), cum), length(cum))
        k99 = something(findfirst(>=(0.99), cum), length(cum))

        results[k] = join([fam_name, n, inst, seed, m, ramp.id,
                           E_cc, k_match, k90, k99, length(cum)], ",")
        println("done: $fam_name n=$n inst=$inst ramp=$(ramp.id)  " *
                "E_cc=$(round(E_cc, digits=4)) k_match=$k_match k90=$k90 k99=$k99")
    end

    outpath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,ramp,E_cc,k_match,k90,k99,rank_max")
        foreach(r -> println(io, r), results)
    end
    println("\nE2.3 sweep complete → $outpath")
end

main()
