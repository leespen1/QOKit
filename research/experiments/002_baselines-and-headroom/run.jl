#=
E1.1 — How much headroom above trivial baselines does proxy parameter setting
actually have, per graph family?

For each instance we compute approximation ratios (AR = ⟨C⟩/c_opt) of:
  (a) uniform random measurement        — ⟨C⟩ = mean(costs)
  (b) random balanced partition         — exact mean over all balanced bitstrings
  (c) true-QAOA grid ceiling            — best AR over a (γ,β) grid (p=1) and a
                                          linear-ramp endpoint grid (p=3)
  (d) proxy-set parameters              — real-QAOA AR at the argmax of the
                                          exact-compression proxy (same-instance
                                          empirical N and empirical P) over the
                                          same grids

Downstream metrics: regret = (c) − (d), value-added = (d) − (b). If value-added
is ≈ 0 within noise on all families, parameter setting is ill-posed at this
scale and the program reframes (Phase-1 gate).

Run from the repo root (full run ~15–60 min with threads):
  JULIA_NUM_THREADS=auto julia --project research/experiments/002_baselines-and-headroom/run.jl
Smoke test (~1 min): E1_SMOKE=1 julia --project ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads
using Statistics: mean

const SMOKE = get(ENV, "E1_SMOKE", "0") == "1"

const SEED = 20260611
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
    # No commas in family names — they are CSV field values
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("WS(k=4;b=0.5)", (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]

"Exact mean cost over all balanced bitstrings (n even)."
function balanced_mean_cost(costs, n)
    @assert iseven(n) "balanced partitions need even n"
    total = 0.0
    count = 0
    for x in 0:(length(costs) - 1)
        if count_ones(x) == n ÷ 2
            total += costs[x + 1]
            count += 1
        end
    end
    return total / count
end

"Best real-QAOA expectation over a list of (γs, βs) schedules (threaded)."
function real_grid_ceiling(costs, n, schedules)
    vals = zeros(length(schedules))
    @threads for k in eachindex(schedules)
        γs, βs = schedules[k]
        vals[k] = qaoa_expectation(costs, n, γs, βs)
    end
    best = argmax(vals)
    return vals[best], schedules[best]
end

"Argmax of the exact-compression proxy over K schedules given as K×p matrices.
Returns (index, γs, βs) of the best schedule."
function proxy_argmax_schedule(N, P_emp, n, γmat, βmat)
    Qs = QAOA_proxy_multi(N, γmat, βmat)
    vals = vec(expectation(Qs[end], P_emp, n))
    best = argmax(vals)
    return best, γmat[best, :], βmat[best, :]
end

function main()
    outpath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,c_opt,ar_uniform,ar_balanced," *
                    "ar_ceiling_p1,ar_proxy_p1,gamma_p1,beta_p1," *
                    "ar_ceiling_p3,ar_proxy_p3,g1,gf,b1,bf")

        for (fam_idx, (fam_name, gen)) in enumerate(FAMILIES), n in NS, inst in 1:INSTANCES
            seed = SEED + 10_000 * fam_idx + 100 * n + inst
            rng = MersenneTwister(seed)
            edges = gen(rng, n)
            m = length(edges)
            @assert m > 0 "empty graph generated (family=$fam_name, seed=$seed)"
            costs = maxcut_costs(n, edges)
            c_opt = maximum(costs)

            ar_uniform = mean(costs) / c_opt
            ar_balanced = balanced_mean_cost(costs, n) / c_opt

            # Exact-compression proxy ingredients (same-instance empirical N and P)
            N = get_homogeneous_distribution_from_costs_direct(costs, m, n)
            counts = zeros(Int, m + 1)
            for c in costs
                counts[Int(c) + 1] += 1
            end
            P_emp = counts ./ (1 << n)

            # ── p = 1: full (γ, β) grid ──
            p1_schedules = vec([([γ], [β]) for γ in P1_γ, β in P1_β])
            ceil_p1, _ = real_grid_ceiling(costs, n, p1_schedules)
            γmat1 = reshape([γ for γ in P1_γ, β in P1_β] |> vec, :, 1)
            βmat1 = reshape([β for γ in P1_γ, β in P1_β] |> vec, :, 1)
            _, γ1, β1 = proxy_argmax_schedule(N, P_emp, n, γmat1, βmat1)
            proxy_p1 = qaoa_expectation(costs, n, γ1, β1)

            # ── p = 3: linear-ramp endpoint grid ──
            ramp_combos = vec([(g1, gf, b1, bf) for g1 in RAMP_γ, gf in RAMP_γ,
                                                    b1 in RAMP_β, bf in RAMP_β])
            ramp_schedules = [linear_ramp(c..., P_RAMP) for c in ramp_combos]
            ceil_p3, _ = real_grid_ceiling(costs, n, ramp_schedules)
            γmat3, βmat3 = linear_ramp_matrix(
                [c[1] for c in ramp_combos], [c[2] for c in ramp_combos],
                [c[3] for c in ramp_combos], [c[4] for c in ramp_combos], P_RAMP)
            best3, γ3, β3 = proxy_argmax_schedule(N, P_emp, n, γmat3, βmat3)
            proxy_p3 = qaoa_expectation(costs, n, γ3, β3)
            g1, gf, b1, bf = ramp_combos[best3]

            row = join([
                fam_name, n, inst, seed, m, c_opt,
                ar_uniform, ar_balanced,
                ceil_p1 / c_opt, proxy_p1 / c_opt, γ1[1], β1[1],
                ceil_p3 / c_opt, proxy_p3 / c_opt, g1, gf, b1, bf,
            ], ",")
            println(io, row)
            flush(io)
            println("done: $fam_name n=$n inst=$inst  m=$m  " *
                    "AR bal/ceil1/prox1/ceil3/prox3 = " *
                    join(round.([ar_balanced, ceil_p1 / c_opt, proxy_p1 / c_opt,
                                 ceil_p3 / c_opt, proxy_p3 / c_opt], digits=4), " / "))
        end
    end
    println("\nE1.1 sweep complete → $outpath")
end

main()
