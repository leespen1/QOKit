#=
E3.1 — Does the Phase-1 picture survive at scale (n = 16, 18)?

Per instance, mirroring the gate experiment (004) so results are directly
comparable across n:
  - true-QAOA grid ceilings at p=1 (40×40) and p=3 linear ramps (8⁴), on GPU
    when available;
  - parameter-setting AR of the exact-compression proxy (exact empirical N via
    GPU homodist) and of the sampled-N proxy (S = 10 per class, the Phase-3
    scalable route) — their agreement validates sampled-N parameter setting;
  - baselines (uniform, exact balanced-partition mean);
  - leakage diagnostics (Σλ, overlap, distance) along the sampled-N-chosen
    p=3 schedule via the exact compressed trajectory.

Slurm array: tasks 1..14 map to (family, n) pairs, 20 instances each
(SLURM_ARRAY_TASK_ID; unset/0 = run all tasks serially, for local smoke).
Each task writes results_task<ID>.csv; concatenate for analysis.

Submit:  sbatch research/experiments/010_scaleup-ranking/run.sb
Smoke:   E1_SMOKE=1 julia --project research/experiments/010_scaleup-ranking/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads
using Statistics: mean

const USE_GPU = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end
println("USE_GPU = ", USE_GPU)

const SMOKE = get(ENV, "E1_SMOKE", "0") == "1"

const SEED = 20260611   # same scheme as experiments 002/004-009
const NS = SMOKE ? [12] : [16, 18]
const INSTANCES = SMOKE ? 2 : 20
const S_N = 10
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
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("WS(k=4;b=0.5)", (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]

const TASKS = [(fam_idx, n) for n in NS for fam_idx in eachindex(FAMILIES)]

"Best real-QAOA expectation over schedules: GPU if available, else threaded CPU."
function grid_ceiling(costs, n, schedules)
    if USE_GPU
        costs_dev = CUDA.CuArray(costs)
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

"Real-QAOA AR at the argmax of a proxy over K schedules (K×p matrices)."
function proxy_choice_ar(N, P, n, γmat, βmat, schedules, costs, c_opt)
    Qs = QAOA_proxy_multi(N, γmat, βmat)
    best = argmax(vec(expectation(Qs[end], P, n)))
    γs, βs = schedules[best]
    ar = if USE_GPU
        gpu_qaoa_expectation_batched(CUDA.CuArray(costs), n, γs, βs) / c_opt
    else
        qaoa_expectation(costs, n, γs, βs) / c_opt
    end
    return ar, γs, βs
end

function balanced_mean_cost(costs, n)
    @assert iseven(n) "balanced partitions need even n"
    total = 0.0; count = 0
    for x in 0:(length(costs) - 1)
        if count_ones(x) == n ÷ 2
            total += costs[x + 1]; count += 1
        end
    end
    return total / count
end

function main()
    tid = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
    tasks = tid == 0 ? TASKS : [TASKS[tid]]
    println("running tasks: ", tasks)

    p1_schedules = vec([([γ], [β]) for γ in P1_γ, β in P1_β])
    γmat1 = reshape([s[1][1] for s in p1_schedules], :, 1)
    βmat1 = reshape([s[2][1] for s in p1_schedules], :, 1)
    ramp_combos = vec([(g1, gf, b1, bf) for g1 in RAMP_γ, gf in RAMP_γ,
                                            b1 in RAMP_β, bf in RAMP_β])
    ramp_schedules = [linear_ramp(c..., P_RAMP) for c in ramp_combos]
    γmat3, βmat3 = linear_ramp_matrix(
        [c[1] for c in ramp_combos], [c[2] for c in ramp_combos],
        [c[3] for c in ramp_combos], [c[4] for c in ramp_combos], P_RAMP)

    rows = String[]
    for (fam_idx, n) in tasks
        fam_name, gen = FAMILIES[fam_idx]
        for inst in 1:INSTANCES
            seed = SEED + 10_000 * fam_idx + 100 * n + inst
            edges = gen(MersenneTwister(seed), n)
            m = length(edges)
            costs = maxcut_costs(n, edges)
            c_opt = maximum(costs)
            ar_uniform = mean(costs) / c_opt
            ar_balanced = balanced_mean_cost(costs, n) / c_opt

            counts = zeros(Int, m + 1)
            for c in costs
                counts[Int(c) + 1] += 1
            end
            P_emp = counts ./ (1 << n)

            N_samp = sampled_homogeneous_distribution(
                costs, m, n; samples_per_class=S_N, rng=MersenneTwister(seed + 777))
            N_exact = if USE_GPU
                Array(gpu_get_homogeneous_distribution_from_costs_direct(costs, m, n))
            elseif n <= 14
                get_homogeneous_distribution_from_costs_direct(costs, m, n)
            else
                nothing   # O(4^n) CPU is out of budget; sampled-only
            end

            row = Any[fam_name, n, inst, seed, m, ar_uniform, ar_balanced]
            chosen3 = nothing
            for (γmat, βmat, schedules, p) in ((γmat1, βmat1, p1_schedules, 1),
                                               (γmat3, βmat3, ramp_schedules, P_RAMP))
                ceiling = grid_ceiling(costs, n, schedules) / c_opt
                ar_samp, γs, βs = proxy_choice_ar(N_samp, P_emp, n, γmat, βmat,
                                                  schedules, costs, c_opt)
                p == P_RAMP && (chosen3 = (γs, βs))
                ar_exact = if N_exact === nothing
                    NaN
                else
                    proxy_choice_ar(N_exact, P_emp, n, γmat, βmat,
                                    schedules, costs, c_opt)[1]
                end
                row = vcat(row, [ceiling, ar_samp, ar_exact])
            end

            traj = compressed_qaoa_trajectory(costs, n, chosen3...; num_costs=m + 1)
            row = vcat(row, [sum(traj.leakage), traj.overlap[end], traj.distance[end]])

            push!(rows, join(row, ","))
            println("done: $fam_name n=$n inst=$inst m=$m")
            flush(stdout)
        end
    end

    suffix = tid == 0 ? (SMOKE ? "_smoke" : "_all") : "_task$tid"
    outpath = joinpath(@__DIR__, "results$suffix.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,ar_uniform,ar_balanced," *
                    "ceil_p1,ar_samp_p1,ar_exact_p1," *
                    "ceil_p3,ar_samp_p3,ar_exact_p3," *
                    "sum_leakage_p3,overlap_p3,distance_p3")
        foreach(r -> println(io, r), rows)
    end
    println("E3.1 task(s) complete → $outpath")
end

main()
