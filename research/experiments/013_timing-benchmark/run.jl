#=
E013 — Timing benchmark: what does each stage of parameter setting cost?

One honest wall-clock table for the paper. Per n, on one G(n,0.5) instance
(fixed seed), median over reps of:

  distributions
    exact_N      — empirical N(v';d,v): GPU direct homodist, O(4^n)
                   (CPU direct for n ≤ 12 when no GPU, e.g. smoke)
    sampled_N    — stratified sampled N, S=10 per class, O(S m 2^n), CPU
  proxy parameter setting (given N; includes expectation + argmax)
    proxy_sweep_p1   — all 40×40 (γ,β) schedules via QAOA_proxy_multi
    proxy_sweep_p3   — all 8^4 linear-ramp endpoint schedules at p=3
    proxy_sweep_p20  — the same 8^4 schedules at p=20
  real-QAOA reference (statevector; GPU when available, costs pre-uploaded)
    sv_single_p1     — one (γ,β) expectation
    sv_single_p20    — one 20-layer ramp expectation
    ceiling_p1       — brute-force grid ceiling, 1600 evaluations
    ceiling_p3       — brute-force ramp ceiling, 4096 evaluations

Submit:  sbatch research/experiments/013_timing-benchmark/run.sb
Smoke:   E13_SMOKE=1 julia --project research/experiments/013_timing-benchmark/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister
using Statistics: median
using Base.Threads: @threads

const USE_GPU = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end
println("USE_GPU = ", USE_GPU)
println("CPU: ", Sys.cpu_info()[1].model, " (", Sys.CPU_THREADS, " threads, ",
        Threads.nthreads(), " Julia threads)")
USE_GPU && println("GPU: ", CUDA.name(CUDA.device()))

const SMOKE = get(ENV, "E13_SMOKE", "0") == "1"
const SEED = 20260612
const NS = SMOKE ? [10, 12] : [12, 14, 16, 18, 20]
const S_N = 10
const P1_GRID_LEN = SMOKE ? 12 : 40
const RAMP_GRID_LEN = SMOKE ? 4 : 8
const P_RAMP = 3
const P_DEEP = 20
const REPS = SMOKE ? 2 : 3

const P1_γ = collect(range(0.0, π; length=P1_GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=P1_GRID_LEN))
const RAMP_γ = collect(range(0.05, 1.6; length=RAMP_GRID_LEN))
const RAMP_β = collect(range(0.05, 0.8; length=RAMP_GRID_LEN))

"Median wall-clock of f() over up to `reps` runs after a warmup call;
stops early once `budget` seconds of measured time accumulate."
function timeit(f; reps=REPS, budget=300.0)
    f()
    ts = Float64[]
    for _ in 1:reps
        t = @elapsed f()
        push!(ts, t)
        sum(ts) > budget && break
    end
    return median(ts), length(ts)
end

"Full proxy parameter-setting decision: sweep all schedules, pick argmax."
function proxy_sweep(N, P_emp, n, γmat, βmat)
    Qs = QAOA_proxy_multi(N, γmat, βmat)
    return argmax(vec(expectation(Qs[end], P_emp, n)))
end

"Brute-force ceiling over schedules (statevector per schedule)."
function ceiling(costs, costs_dev, n, schedules)
    if USE_GPU
        best = -Inf
        for (γs, βs) in schedules
            best = max(best, gpu_qaoa_expectation_batched(costs_dev, n, γs, βs))
        end
        return best
    end
    vals = zeros(length(schedules))
    @threads for k in eachindex(schedules)
        γs, βs = schedules[k]
        vals[k] = qaoa_expectation(costs, n, γs, βs)
    end
    return maximum(vals)
end

sv_single(costs, costs_dev, n, γs, βs) =
    USE_GPU ? gpu_qaoa_expectation_batched(costs_dev, n, γs, βs) :
              qaoa_expectation(costs, n, γs, βs)

function main()
    p1_schedules = vec([([γ], [β]) for γ in P1_γ, β in P1_β])
    γmat1 = reshape([s[1][1] for s in p1_schedules], :, 1)
    βmat1 = reshape([s[2][1] for s in p1_schedules], :, 1)
    ramp_combos = vec([(g1, gf, b1, bf) for g1 in RAMP_γ, gf in RAMP_γ,
                                            b1 in RAMP_β, bf in RAMP_β])
    ramp_schedules = [linear_ramp(c..., P_RAMP) for c in ramp_combos]
    ramp_mats(p) = linear_ramp_matrix(
        [c[1] for c in ramp_combos], [c[2] for c in ramp_combos],
        [c[3] for c in ramp_combos], [c[4] for c in ramp_combos], p)
    γmat3, βmat3 = ramp_mats(P_RAMP)
    γmat20, βmat20 = ramp_mats(P_DEEP)
    γ20, β20 = linear_ramp(0.1, 1.0, 0.5, 0.05, P_DEEP)

    rows = String[]
    rec(n, m, op, t, r) = push!(rows, join((n, m, op, t, r), ","))

    for n in NS
        edges = erdos_renyi_edges(n, 0.5; rng=MersenneTwister(SEED + n))
        m = length(edges)
        costs = maxcut_costs(n, edges)
        costs_dev = USE_GPU ? CUDA.CuArray(costs) : nothing
        counts = zeros(Int, m + 1)
        for c in costs
            counts[Int(c) + 1] += 1
        end
        P_emp = counts ./ (1 << n)

        sampN() = sampled_homogeneous_distribution(
            costs, m, n; samples_per_class=S_N, rng=MersenneTwister(SEED + 777))
        exactN = if USE_GPU
            () -> Array(gpu_get_homogeneous_distribution_from_costs_direct(costs, m, n))
        elseif n <= 12
            () -> get_homogeneous_distribution_from_costs_direct(costs, m, n)
        else
            nothing
        end

        if exactN !== nothing
            t, r = timeit(exactN)
            rec(n, m, "exact_N", t, r)
        end
        t, r = timeit(sampN)
        rec(n, m, "sampled_N", t, r)

        N = exactN === nothing ? sampN() : exactN()
        for (op, γmat, βmat) in (("proxy_sweep_p1", γmat1, βmat1),
                                 ("proxy_sweep_p3", γmat3, βmat3),
                                 ("proxy_sweep_p20", γmat20, βmat20))
            t, r = timeit(() -> proxy_sweep(N, P_emp, n, γmat, βmat))
            rec(n, m, op, t, r)
        end

        t, r = timeit(() -> sv_single(costs, costs_dev, n, [0.4], [0.3]))
        rec(n, m, "sv_single_p1", t, r)
        t, r = timeit(() -> sv_single(costs, costs_dev, n, γ20, β20))
        rec(n, m, "sv_single_p20", t, r)
        t, r = timeit(() -> ceiling(costs, costs_dev, n, p1_schedules))
        rec(n, m, "ceiling_p1", t, r)
        t, r = timeit(() -> ceiling(costs, costs_dev, n, ramp_schedules))
        rec(n, m, "ceiling_p3", t, r)

        println("done: n=$n m=$m")
        flush(stdout)
    end

    suffix = SMOKE ? "_smoke" : ""
    outpath = joinpath(@__DIR__, "results$suffix.csv")
    open(outpath, "w") do io
        println(io, "n,m,op,seconds,reps")
        foreach(r -> println(io, r), rows)
    end
    println("E013 complete → $outpath")
end

main()
