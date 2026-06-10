#=
E1.3 follow-up — Does the norm filter make the analytical PaperProxy safe on
every family, including dense ER(0.5)?

Experiment 004 found that PaperProxy's argmax on dense ER(0.5) lands on
unphysical norm-inflated peaks (predicted ⟨C⟩ exceeding the edge count), while
on other families it is competitive. Theorem 1 implies the exact compression
is contractive, so model-induced norm inflation ‖φ‖² > 1 certifies model error
at zero cost. Here the paper arm is re-run on the full 004 instance set with
grid points filtered to ‖φ‖² ≤ NORM_TOL before taking the argmax.

Run from the repo root (~1 h with threads):
  JULIA_NUM_THREADS=auto julia --project research/experiments/005_norm-filtered-paper-proxy/run.jl
Smoke test: E1_SMOKE=1 julia --project ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E1_SMOKE", "0") == "1"

const SEED = 20260611   # shared with experiments 002/004 → identical instances
const NS = SMOKE ? [10] : [12, 14]
const INSTANCES = SMOKE ? 2 : 30
const P1_GRID_LEN = SMOKE ? 12 : 40
const RAMP_GRID_LEN = SMOKE ? 4 : 8
const P_RAMP = 3
const NORM_TOL = 1.05

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

"Proxy objective values and model-internal state norms over K schedules."
function proxy_vals_norms(N, P, n, γmat, βmat)
    Q = QAOA_proxy_multi(N, γmat, βmat)[end]
    vals = vec(expectation(Q, P, n))
    norm2 = vec(2.0^n .* sum(abs2.(Q) .* P, dims=1))
    return vals, norm2
end

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
            vals, norm2 = proxy_vals_norms(N_paper, P_paper, n, γmat, βmat)

            i_raw = argmax(vals)
            ar_raw = qaoa_expectation(costs, n, schedules[i_raw]...) / c_opt

            # Mild norm drift is ubiquitous for the analytical model (on sparse
            # ramp grids every point can exceed the tolerance); the artifact we
            # must exclude is gross inflation. If no point passes the tolerance,
            # fall back to the least-inflated points (within 10% of min norm²)
            # and record the situation via frac_sane = 0.
            ar_at_filtered_argmax = function (tol)
                sane = norm2 .<= tol
                any(sane) || (sane = norm2 .<= 1.1 * minimum(norm2))
                filtered = copy(vals)
                filtered[.!sane] .= -Inf
                return qaoa_expectation(costs, n, schedules[argmax(filtered)]...) / c_opt
            end
            ar_strict = ar_at_filtered_argmax(NORM_TOL)   # 1.05: any drift rejected
            ar_loose = ar_at_filtered_argmax(2.0)         # only gross inflation rejected
            frac_sane = count(norm2 .<= NORM_TOL) / length(norm2)

            row = vcat(row, [ceiling, ar_raw, ar_strict, ar_loose, frac_sane])
        end
        results[k] = join(row, ",")
        println("done: $(j.fam_name) n=$n inst=$(j.inst)")
    end

    outpath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m," *
                    "ceil_p1,ar_paper_raw_p1,ar_paper_strict_p1,ar_paper_loose_p1,frac_sane_p1," *
                    "ceil_p3,ar_paper_raw_p3,ar_paper_strict_p3,ar_paper_loose_p3,frac_sane_p3")
        foreach(r -> println(io, r), results)
    end
    println("\nE1.3-followup sweep complete → $outpath")
end

main()
