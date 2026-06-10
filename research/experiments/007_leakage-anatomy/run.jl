#=
E2.1 — Leakage anatomy: how does one-layer leakage vary over the (γ, β) plane,
per graph family?

For each grid point we measure two model-independent quantities of the exact
compression:
  - λ_uniform(γ, β): leakage of one QAOA layer applied to the uniform state
    (= λ_1, the first-layer term of Theorem 2);
  - η_F(γ, β) = sqrt(Σ_{c'} λ(c')²): Frobenius operator leakage, where λ(c') is
    the leakage of one layer applied to the normalized cost-class state |c'⟩
    (an upper-bound profile independent of the state).

Questions this feeds:
  - Is the proxy's trust region (small leakage) predictable a priori, and does
    the analytical proxy's ER(0.5) spurious peak (γ≈0.64, β≈1.17, exp 004-006)
    sit in a high-leakage region? If so, leakage itself defines the guard that
    value/norm filters could not provide.
  - Does the density effect (H-density, exp 003) persist when γ is rescaled by
    the edge count (compare grids in raw γ vs γ·m units)?
  - Small-angle scaling: Theorem 3's perturbative prediction λ ~ O(|β|·…).

Run from the repo root (~1 h with threads):
  JULIA_NUM_THREADS=auto julia --project research/experiments/007_leakage-anatomy/run.jl
Smoke test: E1_SMOKE=1 julia --project ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E1_SMOKE", "0") == "1"

const SEED = 20260611   # shared instance seeds with experiments 002/004-006
const NS = SMOKE ? [10] : [12, 14]
const INSTANCES = SMOKE ? 2 : 10
const GRID_LEN = SMOKE ? 8 : 24

const ΓS = collect(range(π / GRID_LEN / 4, π; length=GRID_LEN))
const ΒS = collect(range(π / GRID_LEN / 8, π / 2; length=GRID_LEN))

const FAMILIES = [
    ("ER(0.5)",       (rng, n) -> erdos_renyi_edges(n, 0.5; rng)),
    ("ER(0.25)",      (rng, n) -> erdos_renyi_edges(n, 0.25; rng)),
    ("BA(k=2)",       (rng, n) -> barabasi_albert_edges(n, 2; rng)),
    ("BA(k=4)",       (rng, n) -> barabasi_albert_edges(n, 4; rng)),
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("WS(k=4;b=0.5)", (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]

"One-layer leakage of an arbitrary state (returns residual norm)."
function one_layer_leakage!(work, state, costs, n, γ, β, num_costs)
    copyto!(work, state)
    apply_phase_gate!(work, costs, γ)
    apply_x_mixer!(work, β, n)
    return project_onto_cost_classes(work, costs; num_costs).residual_norm
end

function main()
    jobs = [(fam_idx, fam_name, gen, n, inst)
            for (fam_idx, (fam_name, gen)) in enumerate(FAMILIES)
            for n in NS, inst in 1:INSTANCES]
    results = Vector{String}(undef, length(jobs))

    @threads for k in eachindex(jobs)
        fam_idx, fam_name, gen, n, inst = jobs[k]
        seed = SEED + 10_000 * fam_idx + 100 * n + inst
        edges = gen(MersenneTwister(seed), n)
        m = length(edges)
        costs = maxcut_costs(n, edges)
        num_costs = m + 1

        # Normalized cost-class states |c'⟩ (attained classes only) and |+⟩^n
        counts = zeros(Int, num_costs)
        for c in costs
            counts[Int(c) + 1] += 1
        end
        attained = findall(>(0), counts)
        class_states = [ComplexF64[Int(costs[x]) + 1 == ci ? 1 / sqrt(counts[ci]) : 0.0
                                   for x in eachindex(costs)] for ci in attained]
        uniform = fill(ComplexF64(1 / sqrt(length(costs))), length(costs))
        work = similar(uniform)

        rows = String[]
        for γ in ΓS, β in ΒS
            λu = one_layer_leakage!(work, uniform, costs, n, γ, β, num_costs)
            η2 = 0.0
            for ψc in class_states
                η2 += one_layer_leakage!(work, ψc, costs, n, γ, β, num_costs)^2
            end
            push!(rows, join([fam_name, n, inst, seed, m, γ, β, λu, sqrt(η2)], ","))
        end
        results[k] = join(rows, "\n")
        println("done: $fam_name n=$n inst=$inst (m=$m)")
    end

    outpath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,gamma,beta,lambda_uniform,eta_F")
        foreach(r -> println(io, r), results)
    end
    println("\nE2.1 sweep complete → $outpath")
end

main()
