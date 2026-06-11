#=
E2.4 — the fitted-shape paradox: does a lower entrywise MSE against the
empirical N(c';d,c) mean better parameter setting? (Old-log H2 says no.)

For every instance (same families/seeds as experiments 002/004, first 10
instances), at p=1 on the full 40×40 (γ,β) grid:

  - build the exact empirical N (class-averaged, slices sum to 2^n) and
    empirical P; the empirical-N proxy is the zero-model-error anchor;
  - fit TriangleProxy and NormalProxy to the empirical N by minimizing the
    slice-normalized entrywise MSE with the same smart-random-search and the
    same init/bounds as python/grips (sendai_opt.fit_proxy_to_real,
    run_proxy_study.proxy_init_and_bounds), keeping snapshots at 25%/50% of
    the iteration budget so each instance carries a spectrum of models from
    coarse to converged;
  - add the analytical PaperProxy N (p_eff = 2m/(n(n-1))) as a second anchor;
  - for every model record three kinds of error against the same instance:
      mse_entry      — slice-normalized entrywise MSE (the fit objective),
      eps_raw/_cal   — dynamics-weighted model error ‖Q₁_model − Q₁_emp‖₂,
                       i.e. the one-layer transfer operator difference applied
                       to the uniform initial state, mean over the grid and at
                       the empirical proxy's argmax; _raw uses the model N as
                       produced (historical usage), _cal first rescales every
                       N(c';:,:) slice to sum 2^n (the exact value for the
                       empirical N — a calibration that needs no instance data);
      regret_raw/_cal — grid ceiling AR minus real-QAOA AR at the model's
                       argmax (all models paired with the empirical P).

The paradox is confirmed if regret is uncorrelated (or anticorrelated) with
mse_entry across the model spectrum while it tracks the dynamics-weighted
error; the H2 check is whether regret improves from the init to the fitted
parameters along each fitting trajectory.

Run from the repo root (~15–30 min with threads):
  JULIA_NUM_THREADS=auto julia --project research/experiments/012_fitted-shape-paradox/run.jl
Smoke test: E1_SMOKE=1 julia --project ...
=#

using JuliaQAOA
using Random: MersenneTwister, randn
using Base.Threads: @threads
using Statistics: mean
using LinearAlgebra: norm

const SMOKE = get(ENV, "E1_SMOKE", "0") == "1"

const SEED = 20260611   # shared with experiments 002/004 → identical instances
const NS = SMOKE ? [10] : [12, 14]
const INSTANCES = SMOKE ? 2 : 10
const P1_GRID_LEN = SMOKE ? 12 : 40
const FIT_MAX_ITER = SMOKE ? 60 : 1000

const P1_γ = collect(range(0.0, π; length=P1_GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=P1_GRID_LEN))

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

"Normalize each N(c';:,:) slice to sum 1; all-zero slices stay zero."
function normalize_slices(N::Array{Float64,3})
    out = copy(N)
    for i in axes(out, 1)
        s = sum(@view out[i, :, :])
        s > 0 && (@view(out[i, :, :]) ./= s)
    end
    return out
end

"Rescale each N(c';:,:) slice to sum 2^n (the exact empirical value)."
function calibrate_slices(N::Array{Float64,3}, n::Integer)
    out = copy(N)
    for i in axes(out, 1)
        s = sum(@view out[i, :, :])
        s > 0 && (@view(out[i, :, :]) .*= (2.0^n / s))
    end
    return out
end

"Slice-normalized entrywise MSE between a model homodist and N_emp_norm
(the latter already slice-normalized). Mirrors python/grips
distribution_mean_squared_error(...; normalize=true)."
function slice_mse(N_model::Array{Float64,3}, N_emp_norm::Array{Float64,3})
    (all(isfinite, N_model) && sum(N_model) > 0) || return 1.0e4  # penalize degenerate params
    return mean(abs2, normalize_slices(N_model) .- N_emp_norm)
end

"""
Smart random search, mirroring python/grips/sendai_opt.fit_proxy_to_real:
perturb around the best point, reuse a perturbation while it helps, shrink the
step SDs after `fail_til_shrink` consecutive failures, stop after
`fail_til_end` failures since the last success. Returns the best params, the
best MSE, and snapshots of the best-so-far at the requested iteration counts
(filled with the final best if the search ended earlier).
"""
function smart_random_search(loss, init_params, bounds, rng;
                             max_iter=FIT_MAX_ITER, fail_til_shrink=4,
                             fail_til_end=50, sd_ratio=0.2, checkpoints=Int[])
    lb = [b[1] for b in bounds]
    ub = [b[2] for b in bounds]
    best = clamp.(Float64.(init_params), lb, ub)
    best_mse = loss(best)
    @assert isfinite(best_mse) "loss at the initial parameters must be finite"
    sds = (ub .- lb) .* sd_ratio
    fails = 0
    fails_since_success = 0
    pert = sds .* randn(rng, length(best))
    snaps = Dict{Int,Tuple{Vector{Float64},Float64}}()
    for i in 1:max_iter
        cand = clamp.(best .+ pert, lb, ub)
        mse = loss(cand)
        if mse < best_mse
            best, best_mse = cand, mse
            fails = 0
            fails_since_success = 0
            # keep the same perturbation: it helped
        else
            fails += 1
            fails_since_success += 1
            pert = sds .* randn(rng, length(best))
        end
        if fails >= fail_til_shrink
            sds .*= 2.0 / 3.0
            fails = 0
            pert = sds .* randn(rng, length(best))
        end
        i in checkpoints && (snaps[i] = (copy(best), best_mse))
        fails_since_success >= fail_til_end && break
    end
    for cp in checkpoints
        haskey(snaps, cp) || (snaps[cp] = (copy(best), best_mse))
    end
    return best, best_mse, snaps
end

"Best real-QAOA expectation over schedules (serial — outer loop is threaded)."
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
        seed = instance_seed(fam_idx, n, inst)
        edges = gen(MersenneTwister(seed), n)
        push!(jobs, (; fam_name, n, inst, seed, edges))
    end

    # Analytical PaperProxy N per distinct (n, m), as in experiment 004.
    nm_pairs = unique((j.n, length(j.edges)) for j in jobs)
    println("distinct (n, m) pairs for PaperProxy N: ", length(nm_pairs))
    paper_cache = Dict{Tuple{Int,Int},Array{Float64,3}}()
    cache_lock = ReentrantLock()
    @threads for (n, m) in nm_pairs
        p_eff = 2m / (n * (n - 1))
        N_paper = cpu_compute_homodist(PaperProxy(m, n, p_eff))
        lock(() -> (paper_cache[(n, m)] = N_paper), cache_lock)
    end
    println("PaperProxy cache built")

    p1_schedules = vec([([γ], [β]) for γ in P1_γ, β in P1_β])
    γmat = reshape([s[1][1] for s in p1_schedules], :, 1)
    βmat = reshape([s[2][1] for s in p1_schedules], :, 1)

    # Empirical N arrays serially up front (internally threaded; cannot nest).
    homodists = Vector{Array{Float64,3}}(undef, length(jobs))
    all_costs = Vector{Vector{Float64}}(undef, length(jobs))
    for k in eachindex(jobs)
        j = jobs[k]
        all_costs[k] = maxcut_costs(j.n, j.edges)
        homodists[k] = get_homogeneous_distribution_from_costs_direct(
            all_costs[k], length(j.edges), j.n)
    end
    println("empirical homodists built")

    cp1, cp2 = round(Int, 0.25 * FIT_MAX_ITER), round(Int, 0.5 * FIT_MAX_ITER)
    results = Vector{String}(undef, length(jobs))
    @threads for k in eachindex(jobs)
        j = jobs[k]
        n, m = j.n, length(j.edges)
        costs = all_costs[k]
        c_opt = maximum(costs)

        N_emp = homodists[k]
        # Every attained cost class averages n(x;d,c) over its members, so each
        # nonzero slice must sum to exactly 2^n.
        for i in axes(N_emp, 1)
            s = sum(@view N_emp[i, :, :])
            @assert s == 0 || isapprox(s, 2.0^n; rtol=1e-12) "empirical slice $(i-1) sums to $s, expected 0 or 2^n"
        end
        N_emp_norm = normalize_slices(N_emp)
        counts = zeros(Int, m + 1)
        for cst in costs
            counts[Int(cst) + 1] += 1
        end
        P_emp = counts ./ (1 << n)
        @assert isapprox(sum(P_emp), 1.0; rtol=1e-12)

        ceiling = grid_ceiling(costs, n, p1_schedules) / c_opt

        Q_emp = QAOA_proxy_multi(N_emp, γmat, βmat)[end]
        i_emp = argmax(vec(expectation(Q_emp, P_emp, n)))

        # Fit both shapes; init/bounds mirror run_proxy_study.proxy_init_and_bounds.
        tri_loss(p) = slice_mse(cpu_compute_homodist(OldTriangleProxy(m, n, p...)), N_emp_norm)
        tri_init = [0.0, 0.0, 1.0, 1.0]
        tri_bounds = [(0.0, n^2 / 3), (-10.0, 10.0), (0.005, 2.0), (0.05, 2.0)]
        tri_best, _, tri_snaps = smart_random_search(
            tri_loss, tri_init, tri_bounds, MersenneTwister(j.seed + 777);
            checkpoints=[cp1, cp2])

        norm_loss(p) = slice_mse(cpu_compute_homodist(NormalProxy(m, n, p...)), N_emp_norm)
        norm_init = [m / 2, 1.0, 1.0]
        norm_bounds = [(0.0, Float64(m)), (0.1, max(10.0, Float64(m))), (0.1, max(10.0, Float64(m)))]
        norm_best, _, norm_snaps = smart_random_search(
            norm_loss, norm_init, norm_bounds, MersenneTwister(j.seed + 778);
            checkpoints=[cp1, cp2])

        models = [
            ("emp",       "-",                                 N_emp),
            ("paper",     "-",                                 paper_cache[(n, m)]),
            ("tri_init",  join(tri_init, ";"),                 cpu_compute_homodist(OldTriangleProxy(m, n, tri_init...))),
            ("tri_25",    join(tri_snaps[cp1][1], ";"),        cpu_compute_homodist(OldTriangleProxy(m, n, tri_snaps[cp1][1]...))),
            ("tri_50",    join(tri_snaps[cp2][1], ";"),        cpu_compute_homodist(OldTriangleProxy(m, n, tri_snaps[cp2][1]...))),
            ("tri_fit",   join(tri_best, ";"),                 cpu_compute_homodist(OldTriangleProxy(m, n, tri_best...))),
            ("norm_init", join(norm_init, ";"),                cpu_compute_homodist(NormalProxy(m, n, norm_init...))),
            ("norm_25",   join(norm_snaps[cp1][1], ";"),       cpu_compute_homodist(NormalProxy(m, n, norm_snaps[cp1][1]...))),
            ("norm_50",   join(norm_snaps[cp2][1], ";"),       cpu_compute_homodist(NormalProxy(m, n, norm_snaps[cp2][1]...))),
            ("norm_fit",  join(norm_best, ";"),                cpu_compute_homodist(NormalProxy(m, n, norm_best...))),
        ]

        rows = String[]
        for (name, paramstr, N_model) in models
            @assert all(isfinite, N_model) "model $name produced non-finite N"
            mse_entry = name == "emp" ? 0.0 : slice_mse(N_model, N_emp_norm)

            local eps_stats = Dict{String,Tuple{Float64,Float64}}()
            local ar = Dict{String,Float64}()
            for (tag, N_used) in (("raw", N_model), ("cal", calibrate_slices(N_model, n)))
                Q = QAOA_proxy_multi(N_used, γmat, βmat)[end]
                colerrs = [norm(@view(Q[:, c]) - @view(Q_emp[:, c])) for c in axes(Q, 2)]
                eps_stats[tag] = (mean(colerrs), colerrs[i_emp])
                i_model = argmax(vec(expectation(Q, P_emp, n)))
                γs, βs = p1_schedules[i_model]
                ar[tag] = qaoa_expectation(costs, n, γs, βs) / c_opt
            end
            if name == "emp"
                @assert eps_stats["raw"][1] == 0.0 "empirical model must have zero dynamics error"
            end

            push!(rows, join(Any[j.fam_name, n, j.inst, j.seed, m, name, paramstr,
                                 mse_entry,
                                 eps_stats["raw"][1], eps_stats["raw"][2],
                                 eps_stats["cal"][1], eps_stats["cal"][2],
                                 ceiling,
                                 ar["raw"], ceiling - ar["raw"],
                                 ar["cal"], ceiling - ar["cal"]], ","))
        end
        results[k] = join(rows, "\n")
        println("done: $(j.fam_name) n=$n inst=$(j.inst)")
    end

    outpath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,model,params,mse_entry," *
                    "eps_raw_mean,eps_raw_argmax,eps_cal_mean,eps_cal_argmax," *
                    "ceil_p1,ar_raw,regret_raw,ar_cal,regret_cal")
        foreach(r -> println(io, r), results)
    end
    println("\nE2.4 sweep complete → $outpath")
end

main()
