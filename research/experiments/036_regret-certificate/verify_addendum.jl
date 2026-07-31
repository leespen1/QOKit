#=
E036 addendum check — machine-verify the unnormalized-objective extra term.

The paper's certificate remark states that for the UNNORMALIZED proxy
objective F̂_raw(θ) = ⟨φ|C|φ⟩/c_opt the error acquires one additional term:

  ⟨ψ|C|ψ⟩ − ⟨φ|C|φ⟩ = [⟨ψ|(C−a)|ψ⟩ − ⟨φ|(C−a)|φ⟩] + a(1−‖φ‖²),  a = ⟨C⟩_ψ,

with 1−‖φ_p‖² = Σ_ℓ λ_ℓ² (Theorem 3's norm identity), and hence

  |F − F̂_raw| ≤ [‖ψ−φ‖(σ_ψ + ‖(C−a)φ‖) + |a|(1−‖φ‖²)] / c_opt.

This was audit-verified algebraically (2026-07-31) but, unlike every other
formal statement, had no machine check. This script asserts, on the standard
7 families x n ∈ {10,12} x 3 instances, at p=1 (9 grid points) and p=3
(8 ramps): (A) the norm identity to 1e-10; (B) the decomposition identity to
1e-10; (C) the bound, strictly, at every point.

Run: julia --project research/experiments/036_regret-certificate/verify_addendum.jl (~2 min)
=#

using JuliaQAOA
using Random: MersenneTwister

const SEED = 20260611
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

function checks(costs, n, γs, βs, copt)
    traj = compressed_qaoa_trajectory(costs, n, γs, βs)
    ψ = qaoa_statevector(costs, n, γs, βs)
    Q, counts = traj.Qs[end], traj.counts

    # (A) norm identity: 1 − ‖φ_p‖² = Σ λ²
    φnorm2 = traj.compressed_norm[end]^2
    @assert abs((1 - φnorm2) - sum(abs2, traj.leakage)) < 1e-10

    # state-side moments
    a  = sum(abs2(ψ[i]) * costs[i] for i in eachindex(ψ))
    σψ = sqrt(max(sum(abs2(ψ[i]) * (costs[i] - a)^2 for i in eachindex(ψ)), 0.0))

    # φ-side (unnormalized) moments
    φC  = sum(counts[ci] * abs2(Q[ci]) * (ci - 1) for ci in eachindex(Q))
    φCa2 = sum(counts[ci] * abs2(Q[ci]) * (ci - 1 - a)^2 for ci in eachindex(Q))
    φCa  = sum(counts[ci] * abs2(Q[ci]) * (ci - 1 - a) for ci in eachindex(Q))

    # (B) decomposition identity: ⟨ψCψ⟩−⟨φCφ⟩ = [⟨ψ(C−a)ψ⟩−⟨φ(C−a)φ⟩] + a(1−‖φ‖²)
    lhs = a - φC
    rhs = (0.0 - φCa) + a * (1 - φnorm2)   # ⟨ψ(C−a)ψ⟩ = 0 by choice of a
    @assert abs(lhs - rhs) < 1e-10

    # (C) the bound: |F − F̂_raw| ≤ [‖ψ−φ‖(σψ + ‖(C−a)φ‖) + |a|(1−‖φ‖²)]/c_opt
    distφ = traj.distance[end]
    lhsF  = abs(a - φC) / copt
    rhsF  = (distφ * (σψ + sqrt(max(φCa2, 0.0))) + abs(a) * (1 - φnorm2)) / copt
    @assert lhsF <= rhsF + 1e-12
    return lhsF, rhsF
end

function main()
    npts = 0; worst = 0.0
    for (fi, (fam, gen)) in enumerate(FAMILIES), n in (10, 12), inst in 1:3
        edges = gen(MersenneTwister(instance_seed(fi, n, inst)), n)
        costs = maxcut_costs(n, edges); copt = maximum(costs)
        for g in range(0.1, 1.2; length=3), b in range(0.1, 0.6; length=3)
            checks(costs, n, [g], [b], copt); npts += 1
        end
        for (g1, gf, b1, bf) in ((0.1,0.8,0.5,0.1), (0.2,1.2,0.6,0.2),
                                 (0.05,0.5,0.3,0.05), (0.3,1.0,0.4,0.1),
                                 (0.1,1.5,0.7,0.05), (0.5,0.5,0.2,0.2),
                                 (0.05,1.0,0.5,0.2), (0.2,0.6,0.6,0.3))
            γs, βs = linear_ramp(g1, gf, b1, bf, 3)
            l, r = checks(costs, n, γs, βs, copt)
            worst = max(worst, l / max(r, 1e-300)); npts += 1
        end
    end
    println("all assertions passed at $npts (instance, schedule) points")
    println("worst p=3 bound utilization |F−F̂_raw|/bound = $(round(worst, digits=4))")
end

main()
