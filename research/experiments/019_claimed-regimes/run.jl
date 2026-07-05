#=
E019 — Do the proxy's claimed regimes actually exist?

E013 (timing) and E018 (transfer calibration) leave the proxy dominated at
every size and depth we tested; the paper's remaining practical claims live
in exactly the regimes not yet tested. This experiment tests them:

Part A — DEPTH on the density spread (tasks 1–6): (family, n) ∈
  {ER(0.5), ER(0.25), 3-regular} × {14, 16}, 10 instances, p ∈ {10, 20}
  linear ramps over the standard 8⁴ endpoint grid. E018's only crack in
  transfer's dominance was dense ER at p=3 growing with n; here depth is
  pushed to where compression error is large and transfer might break.

Part B — SIZE beyond the exact-N wall (tasks 7–10): (family, n) ∈
  {ER(0.5), 3-regular} × {22, 26}, 5 instances, p=1 (40×40 grid) and p=3
  ramps. Exact N is O(4^n)-impossible here — the regime the paper assigns
  to the analytical/sampled models. Ceilings are still computable by GPU
  statevector sweep, so every method gets an honest regret.

Methods per instance (all choosing from the same schedule grid as the
ceiling, so regret ≥ 0 by construction):
  ceiling     — best true ⟨C⟩ over the grid (GPU batched statevector)
  exactN      — exact-compression proxy argmax (Part A only; GPU homodist)
  sampledN    — S=10 stratified sampled N + empirical P
  paper_binP  — analytical PaperProxy (effective edge probability), its
                native binomial P
  paper_empP  — same N, empirical P (continuity with exps 010/017)
  transfer_k  — angles optimized on the true landscape of the k-th n=10
                same-family source (k = 1..3; E018 source-seed convention)
  universal   — one family-agnostic schedule: argmax of the mean normalized
                true ⟨C⟩ over all 7 families × 3 sources at n=10

Output: results_task<ID>.csv, long format (one row per instance × p × method).
Submit:  cd research/experiments/019_claimed-regimes && sbatch run.sb
Smoke:   E19_SMOKE=1 julia --project research/experiments/019_claimed-regimes/run.jl
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

const SMOKE = get(ENV, "E19_SMOKE", "0") == "1"
SMOKE || USE_GPU || error("E019 requires a functional GPU (smoke mode runs on CPU)")

const SEED = 20260705
const SOURCE_SEED = 20260704   # E018's convention, for cross-experiment comparability
const SOURCE_N = 10
const N_SOURCES = SMOKE ? 1 : 3
const S_N = 10

const GRID_LEN = SMOKE ? 10 : 40
const RAMP_LEN = SMOKE ? 4 : 8
const P1_γ = collect(range(0.0, π; length=GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=GRID_LEN))
const RAMP_γ = collect(range(0.05, 1.6; length=RAMP_LEN))
const RAMP_β = collect(range(0.05, 0.8; length=RAMP_LEN))

const ALL_FAMILIES = [
    ("ER(0.5)",       (rng, n) -> erdos_renyi_edges(n, 0.5; rng)),
    ("ER(0.25)",      (rng, n) -> erdos_renyi_edges(n, 0.25; rng)),
    ("BA(k=2)",       (rng, n) -> barabasi_albert_edges(n, 2; rng)),
    ("BA(k=4)",       (rng, n) -> barabasi_albert_edges(n, 4; rng)),
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("WS(k=4;b=0.5)", (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]
gen_for(fam) = ALL_FAMILIES[findfirst(t -> t[1] == fam, ALL_FAMILIES)][2]
famindex(fam) = findfirst(==(fam), sort([t[1] for t in ALL_FAMILIES]))

# task table: (part, family, n, instances, depths)
const TASKS = [
    (:A, "ER(0.5)",   14, 10, [10, 20]),
    (:A, "ER(0.25)",  14, 10, [10, 20]),
    (:A, "3-regular", 14, 10, [10, 20]),
    (:A, "ER(0.5)",   16, 10, [10, 20]),
    (:A, "ER(0.25)",  16, 10, [10, 20]),
    (:A, "3-regular", 16, 10, [10, 20]),
    (:B, "ER(0.5)",   22, 5, [1, 3]),
    (:B, "3-regular", 22, 5, [1, 3]),
    (:B, "ER(0.5)",   26, 5, [1, 3]),
    (:B, "3-regular", 26, 5, [1, 3]),
]
const SMOKE_TASK = (:A, "ER(0.5)", 12, 2, [1, 5])

"maxcut_costs, threaded (n=26 is ~10^10 edge checks single-threaded)."
function costs_threaded(n, edges)
    costs = zeros(Float64, 1 << n)
    @threads for x in 0:((1 << n) - 1)
        c = 0
        for (i, j) in edges
            c += ((x >> i) & 1) != ((x >> j) & 1)
        end
        costs[x + 1] = c
    end
    costs
end

"Schedules for depth p: p=1 uses the point grid, else the ramp-endpoint grid."
function schedule_set(p)
    if p == 1
        combos = vec([(γ, β) for γ in P1_γ, β in P1_β])
        scheds = [([c[1]], [c[2]]) for c in combos]
        γmat = reshape([c[1] for c in combos], :, 1)
        βmat = reshape([c[2] for c in combos], :, 1)
        return scheds, γmat, βmat
    end
    combos = vec([(g1, gf, b1, bf) for g1 in RAMP_γ, gf in RAMP_γ,
                                       b1 in RAMP_β, bf in RAMP_β])
    scheds = [linear_ramp(c..., p) for c in combos]
    γmat, βmat = linear_ramp_matrix([c[1] for c in combos], [c[2] for c in combos],
                                    [c[3] for c in combos], [c[4] for c in combos], p)
    return scheds, γmat, βmat
end

"True ⟨C⟩ at one schedule (GPU when available)."
function true_exp(costs, costs_dev, n, γs, βs)
    if USE_GPU
        return gpu_qaoa_expectation_batched(costs_dev, n, γs, βs)
    end
    return qaoa_expectation(costs, n, γs, βs)
end

"Ceiling and per-schedule true values are not needed — only the max."
function grid_ceiling(costs, costs_dev, n, scheds)
    best = -Inf
    for (k, (γs, βs)) in enumerate(scheds)
        best = max(best, true_exp(costs, costs_dev, n, γs, βs))
        if USE_GPU && k % 256 == 0
            CUDA.reclaim()
        end
    end
    best
end

"Argmax of a source instance's TRUE landscape over the schedule set (CPU, n=10)."
function source_argmax(costs, n, scheds)
    vals = zeros(length(scheds))
    @threads for k in eachindex(scheds)
        γs, βs = scheds[k]
        vals[k] = qaoa_expectation(costs, n, γs, βs)
    end
    argmax(vals)
end

"Universal schedule: argmax of mean normalized true ⟨C⟩ over all sources."
function universal_argmax(sources_costs, scheds)
    score = zeros(length(scheds))
    for costs in sources_costs
        c_opt = maximum(costs)
        vals = zeros(length(scheds))
        @threads for k in eachindex(scheds)
            γs, βs = scheds[k]
            vals[k] = qaoa_expectation(costs, SOURCE_N, γs, βs) / c_opt
        end
        score .+= vals
    end
    argmax(score)
end

function main()
    tid = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
    part, fam, n, ninst, depths = SMOKE ? SMOKE_TASK : TASKS[tid == 0 ? 1 : tid]
    println("task $tid: part=$part $fam n=$n instances=$ninst depths=$depths")

    # sources: this family's transfer sources + all-family pool for `universal`
    src_costs = Dict{String, Vector{Vector{Float64}}}()
    for (f, gen) in ALL_FAMILIES
        reps = Vector{Vector{Float64}}()
        for rep in 1:N_SOURCES
            seed = SOURCE_SEED + 1000 * famindex(f) + rep
            push!(reps, maxcut_costs(SOURCE_N, gen(MersenneTwister(seed), SOURCE_N)))
        end
        src_costs[f] = reps
    end
    pool = vcat(values(src_costs)...)

    rows = String[]
    emit(inst, seed, m, p, method, ar, ceil) =
        push!(rows, join(Any[fam, n, inst, seed, m, p, method, ar, ceil, ceil - ar], ","))

    for p in depths
        scheds, γmat, βmat = schedule_set(p)
        println("depth p=$p: $(length(scheds)) schedules")
        transfer_idx = [source_argmax(src_costs[fam][r], SOURCE_N, scheds)
                        for r in 1:N_SOURCES]
        universal_idx = universal_argmax(pool, scheds)
        flush(stdout)

        for inst in 1:ninst
            seed = SEED + 10_000 * famindex(fam) + 100 * n + inst
            edges = gen_for(fam)(MersenneTwister(seed), n)
            m = length(edges)
            costs = costs_threaded(n, edges)
            c_opt = maximum(costs)
            costs_dev = USE_GPU ? CUDA.CuArray(costs) : nothing

            ceil_ar = grid_ceiling(costs, costs_dev, n, scheds) / c_opt
            emit(inst, seed, m, p, "ceiling", ceil_ar, ceil_ar)

            ar_at(k) = true_exp(costs, costs_dev, n, scheds[k]...) / c_opt
            checked(name, k) = begin
                ar = ar_at(k)
                @assert ceil_ar - ar > -1e-8 "$name beat the ceiling"
                emit(inst, seed, m, p, name, ar, ceil_ar)
            end

            # distribution-driven methods
            counts = zeros(Int, m + 1)
            for c in costs
                counts[Int(c) + 1] += 1
            end
            P_emp = counts ./ (1 << n)

            Ns = Vector{Tuple{String, Array{Float64, 3}, Vector{Float64}}}()
            if part == :A
                N_exact = SMOKE ? get_homogeneous_distribution_from_costs_direct(costs, m, n) :
                                  Array(gpu_get_homogeneous_distribution_from_costs_direct(costs, m, n))
                push!(Ns, ("exactN", N_exact, P_emp))
            end
            N_samp = sampled_homogeneous_distribution(costs, m, n;
                         samples_per_class=S_N, rng=MersenneTwister(seed + 777))
            push!(Ns, ("sampledN", N_samp, P_emp))
            paper = PaperProxy(m, n, 2m / (n * (n - 1)))
            N_paper = cpu_compute_homodist(paper)
            P_bin = [P_cost_distribution(paper, c) for c in 0:m]
            push!(Ns, ("paper_binP", N_paper, P_bin))
            push!(Ns, ("paper_empP", N_paper, P_emp))

            for (name, N, P) in Ns
                Qs = QAOA_proxy_multi(N, γmat, βmat)
                checked(name, argmax(vec(expectation(Qs[end], P, n))))
            end
            for r in 1:N_SOURCES
                checked("transfer_$r", transfer_idx[r])
            end
            checked("universal", universal_idx)

            println("done: inst=$inst m=$m p=$p ceil=$(round(ceil_ar; digits=4))")
            flush(stdout)
            USE_GPU && (costs_dev = nothing; CUDA.reclaim())
        end
    end

    suffix = SMOKE ? "_smoke" : "_task$tid"
    outpath = joinpath(@__DIR__, "results$suffix.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,p,method,ar,ceil,regret")
        foreach(r -> println(io, r), rows)
    end
    println("E019 task complete → $outpath")
end

main()
