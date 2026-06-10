#=
E0.1 — Is the homogeneous proxy numerically identical to the compressed
statevector evolution (Theorem 1)?

For ER(0.5) MaxCut instances, compare two independently-computed trajectories:
  (a) QAOA_proxy_basic / QAOA_proxy_single run with the same-instance empirical
      homogeneous distribution N(c'; d, c);
  (b) the true statevector evolved one layer at a time, orthogonally projected
      onto cost classes after every layer (compressed_qaoa_trajectory).
Theorem 1 says these agree exactly on attained cost classes. Any deviation
beyond floating-point noise (1e-10) falsifies the formalization or reveals a
convention mismatch.

Also checks the pi_units convention used by the figure scripts: angles passed
as multiples of π with pi_units=true must reproduce the raw-radians result.

Run from the repo root:  julia --project research/experiments/001_proxy-is-compression/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister

const SEED = 20260610
const NS = [8, 10, 12]
const INSTANCES_PER_N = 5
const P = 3

function er_edges(rng, n, prob)
    edges = Tuple{Int,Int}[]
    for i in 0:(n-2), j in (i+1):(n-1)
        rand(rng) < prob && push!(edges, (i, j))
    end
    return edges
end

rng = MersenneTwister(SEED)
worst = Dict("basic" => 0.0, "single" => 0.0, "pi_units" => 0.0)

for n in NS, inst in 1:INSTANCES_PER_N
    edges = er_edges(rng, n, 0.5)
    m = length(edges)
    costs = maxcut_costs(n, edges)

    # Several (γ, β) regimes: small angles (the proxy's home turf), moderate,
    # and adversarially large
    γs = [0.1, 0.6, 2.3][1:P]
    βs = [0.05, 0.4, 1.9][1:P]

    traj = compressed_qaoa_trajectory(costs, n, γs, βs; num_costs=m + 1)
    N = get_homogeneous_distribution_from_costs_direct(costs, m, n)
    Qs_basic = QAOA_proxy_basic(N, γs, βs)
    Qs_single = QAOA_proxy_single(N, γs, βs)
    Qs_pi = QAOA_proxy_basic(N, γs ./ π, βs ./ π; pi_units=true)

    attained = traj.counts .> 0
    for ℓ in 2:(P + 1)   # layer 0 differs on unattained classes by convention
        worst["basic"] = max(worst["basic"], maximum(abs.(traj.Qs[ℓ] .- Qs_basic[ℓ])))
        worst["single"] = max(worst["single"], maximum(abs.(traj.Qs[ℓ] .- Qs_single[ℓ])))
        worst["pi_units"] = max(worst["pi_units"], maximum(abs.(Qs_basic[ℓ] .- Qs_pi[ℓ])))
    end
    # Layer 0 must agree on attained classes
    worst["basic"] = max(worst["basic"], maximum(abs.(traj.Qs[1][attained] .- Qs_basic[1][attained])))

    println("n=$n instance=$inst m=$m  running max deviations: ",
            join(("$k=$(round(v, sigdigits=3))" for (k, v) in sort(collect(worst))), "  "))
end

println()
for (k, v) in sort(collect(worst))
    println(rpad(k, 10), " max |Δ| = ", v)
    if v > 1e-10
        error("E0.1 FAILED: $k deviates by $v > 1e-10 — Theorem-1 formalization or a convention is wrong; stop and investigate.")
    end
end
println("\nE0.1 PASSED: proxy with same-instance empirical N ≡ project∘apply, to 1e-10.")
