#=
E014 (was E2.4) — The fitted-shape paradox: does Theorem 3's f_d(β)-weighted
model error explain why lower entrywise-MSE fits can set parameters WORSE?

Theorem 3 says the proxy update sees the model N(c';d,c) only through the
one-layer transfer contraction
    T(c',c; β) = Σ_d f_d(β) N(c';d,c),   f_d(β) = cos(β)^(n-d) (-i sin β)^d,
so the parameter-setting-relevant model error is the transfer-matrix error
    E_w(β)² = Σ_{c',c} | Σ_d f_d(β) ΔN(c';d,c) |²
(the γ phase factors are unimodular per column and drop out of the Frobenius
norm), NOT the entrywise error ‖ΔN‖_F that the G-RIPS shape fits minimized.

Three parts, all at p=1 on the standard 40×40 (γ,β) grid with true-QAOA grid
ceilings (mirroring experiments 002/010):

A. Controlled perturbations of the exact empirical N at MATCHED entrywise
   error ε‖N‖_F but very different weighted error:
     lowd   — d-profile aligned with Re f(β=0.3) (maximally visible),
     highd  — mass at d ≈ n/2 where f_d(β) is negligible for small β,
     nulled — low-d profile with span{Re f, Im f : β ∈ PERT_β_REF} projected
              out (invisible at those βs by construction),
     random — a random d-profile.
   Prediction: at fixed ε (fixed entrywise MSE), regret tracks E_w and spans
   orders of magnitude across shapes.

B. Actual shape fits as G-RIPS did: TriangleProxy and NormalProxy fit to the
   exact empirical N by entrywise error (random search), plus the analytical
   PaperProxy with effective edge probability. Question: where entrywise MSE
   misranks models against regret, does E_w rank them correctly?

C. The prescriptive fix: fit the same shapes in the weighted norm (E_w over
   reference βs). Prediction: weighted-norm fits set parameters at least as
   well as entrywise fits despite (much) larger entrywise MSE.

Global scale is a gauge freedom of the proxy argmax (a uniform s^{2p} factor
on |Q|² across the parameter grid), so every error metric and every fit
objective first optimizes a global scale s in closed form.

Slurm array: tasks 1..10 map to (family, n) pairs, 15 instances each
(SLURM_ARRAY_TASK_ID; unset/0 = run all tasks serially).
Each task writes results_task<ID>.csv (long format: one row per variant).

Submit:  cd research/experiments/014_fitted-shape-paradox && sbatch run.sb
Smoke:   E14_SMOKE=1 julia --project research/experiments/014_fitted-shape-paradox/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister, randn
using Base.Threads: @threads
using Statistics: mean
using LinearAlgebra: norm, dot, normalize, qr

const USE_GPU = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end
println("USE_GPU = ", USE_GPU)

const SMOKE = get(ENV, "E14_SMOKE", "0") == "1"

const SEED = 20260703
const NS = SMOKE ? [10] : [12, 14]
const INSTANCES = SMOKE ? 2 : 15
const GRID_LEN = SMOKE ? 10 : 40
const N_INIT = SMOKE ? 60 : 1500      # random-search initial samples
const N_ROUND = SMOKE ? 30 : 400      # samples per refinement round
const ROUNDS = SMOKE ? 2 : 4

const γ_GRID = collect(range(0.0, π; length=GRID_LEN))
const β_GRID = collect(range(0.0, π/2; length=GRID_LEN))
const β_REF_FIT = collect(range(0.05, π/2 - 0.05; length=8))  # weighted-fit βs
const PERT_β_REF = [0.15, 0.3, 0.45]  # null-projection reference βs
const EPSILONS = [0.01, 0.05, 0.2, 0.5]

const FAMILIES = [
    ("ER(0.5)",       (rng, n) -> erdos_renyi_edges(n, 0.5; rng)),
    ("ER(0.25)",      (rng, n) -> erdos_renyi_edges(n, 0.25; rng)),
    ("BA(k=4)",       (rng, n) -> barabasi_albert_edges(n, 4; rng)),
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]

const TASKS = [(fam_idx, n) for n in NS for fam_idx in eachindex(FAMILIES)]

# ---------------------------------------------------------------- metrics ---

"Stack N over d: (n+1) × (m+1)² matrix D with D[1+d, j] = N[c'(j), 1+d, c(j)]."
flatten_dslices(N) = reshape(permutedims(N, (2, 1, 3)), size(N, 2), :)

"β-factor matrix F, (n+1) × K, with F[1+d, k] = f_d(β_k)."
βfactor_matrix(βs, n) = get_β_factors(collect(βs), n)

"Weighted transfer rows: K × (m+1)² complex, row k = vec of Σ_d f_d(β_k) N."
weighted_rows(N, F) = transpose(F) * flatten_dslices(N)

"Gauge-fixed relative entrywise error: (min_s ‖sA−B‖_F / ‖B‖_F, s*)."
function gauge_relerr(A, B)
    aa = dot(vec(A), vec(A))
    aa == 0 && return 1.0, 0.0
    s = dot(vec(A), vec(B)) / aa
    return norm(s .* A .- B) / norm(B), s
end

"Gauge-fixed relative weighted error over the βs in F: (relerr, s*)."
function gauge_weighted_relerr(A, B, F)
    TA = weighted_rows(A, F)
    TB = weighted_rows(B, F)
    aa = real(dot(TA, TA))
    aa == 0 && return 1.0, 0.0
    s = real(dot(TA, TB)) / aa   # best real global scale
    return norm(s .* TA .- TB) / norm(TB), s
end

# --------------------------------------------------------- proxy machinery ---

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

#=
Two rankings of the same proxy sweep:
  raw — Eq. 9 of the paper, ⟨C⟩ = Σ_c 2^n P(c)|Q(c)|² c, the convention of
        experiments 001–010. Vulnerable to norm inflation: model error can
        pump ‖Q‖ at specific (γ,β) and hijack the argmax (the exp-004–006
        pathology), which smoke runs showed happens for ANY coherent ΔN.
  nrm — the same quantity divided by the state weight Σ_c 2^n P(c)|Q(c)|²,
        i.e. ⟨C⟩ of the normalized compressed state. Gauge-invariant per
        parameter set; immune to pure norm inflation. This is the instrument
        that isolates how model error moves the LANDSCAPE SHAPE, and doubles
        as a test of "normalize, don't veto" (exps 005/006 showed value
        vetoes fail; division may be the recipe that works).
=#
function proxy_argmaxes(N, P, n, γmat, βmat)
    Qs = QAOA_proxy_multi(N, γmat, βmat)
    vals = vec(expectation(Qs[end], P, n))
    weights = vec(sum((2.0^n .* P) .* abs2.(Qs[end]); dims=1))
    best_raw = argmax(vals)
    nrm = [w > 1e-300 ? v / w : -Inf for (v, w) in zip(vals, weights)]
    best_nrm = argmax(nrm)
    return best_raw, vals[best_raw], weights[best_raw], best_nrm
end

function real_ar(costs, n, γs, βs, c_opt)
    if USE_GPU
        return gpu_qaoa_expectation_batched(CUDA.CuArray(costs), n, γs, βs) / c_opt
    end
    return qaoa_expectation(costs, n, γs, βs) / c_opt
end

# ----------------------------------------------------------- perturbations ---

"Four unit d-profiles: aligned low-d, high-d, nulled-at-reference-βs, random."
function perturbation_profiles(n, rng)
    f_mid = βfactor_matrix([0.3], n)[:, 1]
    lowd = normalize(real.(f_mid))

    highd = zeros(n + 1)
    highd[(n ÷ 2):(n ÷ 2 + 2)] .= 1 / sqrt(3)   # d = n/2−1, n/2, n/2+1

    F_ref = βfactor_matrix(PERT_β_REF, n)
    basis = Matrix{Float64}(undef, n + 1, 0)
    for k in axes(F_ref, 2)
        basis = hcat(basis, real.(F_ref[:, k]), imag.(F_ref[:, k]))
    end
    Q = Matrix(qr(basis).Q)[:, 1:size(basis, 2)]
    # start from a low-d bump NOT in the reference span (lowd itself is in it)
    bump = normalize(exp.(-(0:n) ./ 2))
    v = bump - Q * (Q' * bump)
    @assert norm(v) > 1e-10 "null projection annihilated the low-d bump"
    nulled = normalize(v)

    return [("lowd", lowd), ("highd", highd),
            ("nulled", nulled), ("random", normalize(randn(rng, n + 1)))]
end

"ΔN[c',d,c] = u[1+d] · ‖N[c',:,c]‖, rescaled so ‖ΔN‖_F = ε ‖N‖_F."
function build_perturbation(N, u, ε)
    m1, n1, _ = size(N)
    w = [norm(@view N[i, :, j]) for i in 1:m1, j in 1:m1]
    ΔN = [u[d] * w[i, j] for i in 1:m1, d in 1:n1, j in 1:m1]
    ΔN .*= ε * norm(N) / norm(ΔN)
    return ΔN
end

# ------------------------------------------------------------------ fitting ---

struct ParamSpec
    lo::Float64
    hi::Float64
    logscale::Bool
end

function sample_uniform(rng, s::ParamSpec)
    s.logscale && return exp(log(s.lo) + rand(rng) * (log(s.hi) - log(s.lo)))
    return s.lo + rand(rng) * (s.hi - s.lo)
end

"Sample in a shrunken box around x0 (log-space box for logscale params)."
function sample_near(rng, s::ParamSpec, x0, shrink)
    if s.logscale
        w = shrink * (log(s.hi) - log(s.lo))
        return exp(clamp(log(x0) + (rand(rng) - 0.5) * w, log(s.lo), log(s.hi)))
    end
    w = shrink * (s.hi - s.lo)
    return clamp(x0 + (rand(rng) - 0.5) * w, s.lo, s.hi)
end

"""
Random search: N_INIT uniform samples, then ROUNDS shrinking-box refinements
of N_ROUND samples each. Candidates are drawn serially from rng (deterministic
for any thread count); objective evaluations are threaded. build_N(x) may
throw (e.g. non-PD covariance) — such candidates score Inf.
"""
function random_search(build_N, specs, objective, rng)
    score = function (x)
        Nfit = try
            build_N(x)
        catch
            return Inf
        end
        all(isfinite, Nfit) || return Inf
        return objective(Nfit)
    end
    best_x = [sample_uniform(rng, s) for s in specs]
    best_v = score(best_x)
    function consume!(X)
        vals = fill(Inf, length(X))
        @threads for i in eachindex(X)
            vals[i] = score(X[i])
        end
        for i in eachindex(X)
            if vals[i] < best_v
                best_v = vals[i]
                best_x = X[i]
            end
        end
    end
    consume!([[sample_uniform(rng, s) for s in specs] for _ in 1:N_INIT])
    for r in 1:ROUNDS
        shrink = 0.5^r
        consume!([[sample_near(rng, specs[k], best_x[k], shrink)
                   for k in eachindex(specs)] for _ in 1:N_ROUND])
    end
    @assert isfinite(best_v) "random search found no finite-objective candidate"
    return best_x, best_v
end

triangle_specs() = [ParamSpec(0.05, 8.0, true),    # height_adjustment
                    ParamSpec(0.0, 1.0, false),    # center_adjustment
                    ParamSpec(0.02, 0.48, false),  # left_angle
                    ParamSpec(0.02, 0.48, false)]  # right_angle

normal_specs(m) = [ParamSpec(0.2m, 0.8m, false),   # cost_mean
                   ParamSpec(0.01, 10.0m, true),   # cov_1
                   ParamSpec(0.01, 10.0m, true)]   # cov_2

# --------------------------------------------------------------------- main ---

function main()
    tid = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
    tasks = tid == 0 ? TASKS : [TASKS[tid]]
    println("running tasks: ", tasks)

    p1_schedules = vec([([γ], [β]) for γ in γ_GRID, β in β_GRID])
    γmat = reshape([s[1][1] for s in p1_schedules], :, 1)
    βmat = reshape([s[2][1] for s in p1_schedules], :, 1)

    rows = String[]
    for (fam_idx, n) in tasks
        fam_name, gen = FAMILIES[fam_idx]
        F_grid = βfactor_matrix(β_GRID, n)
        F_fit = βfactor_matrix(β_REF_FIT, n)
        for inst in 1:INSTANCES
            seed = SEED + 10_000 * fam_idx + 100 * n + inst
            rng = MersenneTwister(seed)
            edges = gen(rng, n)
            m = length(edges)
            costs = maxcut_costs(n, edges)
            c_opt = maximum(costs)

            counts = zeros(Int, m + 1)
            for c in costs
                counts[Int(c) + 1] += 1
            end
            P_emp = counts ./ (1 << n)

            N_exact = get_homogeneous_distribution_from_costs_direct(costs, m, n)
            for cp in 0:m   # every attained class must see all 2^n bitstrings
                counts[cp + 1] == 0 && continue
                @assert isapprox(sum(@view N_exact[cp + 1, :, :]), 2.0^n; rtol=1e-9)
            end

            ceiling_ar = grid_ceiling(costs, n, p1_schedules) / c_opt

            # reference: exact-N proxy choices under both rankings
            ref_raw, _, _, ref_nrm = proxy_argmaxes(N_exact, P_emp, n, γmat, βmat)
            point(k) = (p1_schedules[k][1][1], p1_schedules[k][2][1])
            refγ_raw, refβ_raw = point(ref_raw)
            refγ_nrm, refβ_nrm = point(ref_nrm)
            F_star_raw = βfactor_matrix([refβ_raw], n)
            F_star_nrm = βfactor_matrix([refβ_nrm], n)

            ar_cache = Dict{Int, Float64}()
            ar_at! = function (k)
                get!(ar_cache, k) do
                    γs, βs = p1_schedules[k]
                    real_ar(costs, n, γs, βs, c_opt)
                end
            end

            emit = function (variant, Nvar, P; fitobj=NaN, fitparams="")
                braw, pred, normw, bnrm = proxy_argmaxes(Nvar, P, n, γmat, βmat)
                ar_raw = ar_at!(braw)
                ar_nrm = ar_at!(bnrm)
                @assert ceiling_ar - ar_raw > -1e-8 "grid point beat the ceiling"
                @assert ceiling_ar - ar_nrm > -1e-8 "grid point beat the ceiling"
                mse_rel, _ = gauge_relerr(Nvar, N_exact)
                ew_avg, _ = gauge_weighted_relerr(Nvar, N_exact, F_grid)
                ew_star_raw, _ = gauge_weighted_relerr(Nvar, N_exact, F_star_raw)
                ew_star_nrm, _ = gauge_weighted_relerr(Nvar, N_exact, F_star_nrm)
                γr, βr = point(braw)
                γn, βn = point(bnrm)
                push!(rows, join(Any[fam_name, n, inst, seed, m, c_opt, ceiling_ar,
                                     variant, mse_rel, ew_avg, ew_star_raw, ew_star_nrm,
                                     ar_raw, ceiling_ar - ar_raw, γr - refγ_raw, βr - refβ_raw,
                                     ar_nrm, ceiling_ar - ar_nrm, γn - refγ_nrm, βn - refβ_nrm,
                                     pred, normw, fitobj, fitparams], ","))
            end

            emit("exactN", N_exact, P_emp)

            # Part A: matched-entrywise-MSE perturbations
            for (shape, u) in perturbation_profiles(n, rng), ε in EPSILONS
                emit("pert:$shape:$ε", N_exact .+ build_perturbation(N_exact, u, ε), P_emp)
            end

            # Part B/C: shape fits in both norms + analytical PaperProxy
            mse_obj = Nfit -> gauge_relerr(Nfit, N_exact)[1]
            wgt_obj = Nfit -> gauge_weighted_relerr(Nfit, N_exact, F_fit)[1]
            builders = [
                ("triangle", x -> cpu_compute_homodist(IntuitiveTriangleProxy(m, n, x...)), triangle_specs()),
                ("normal",   x -> cpu_compute_homodist(NormalProxy(m, n, x...)),            normal_specs(m)),
            ]
            for (model, build, specs) in builders
                for (objname, obj) in (("mse", mse_obj), ("wgt", wgt_obj))
                    x, v = random_search(build, specs, obj,
                                         MersenneTwister(seed + hash(model * objname) % 10_000))
                    emit("fit:$model:$objname", build(x), P_emp;
                         fitobj=v, fitparams=join(round.(x; sigdigits=6), ";"))
                end
            end

            paper = PaperProxy(m, n, 2m / (n * (n - 1)))
            N_paper = cpu_compute_homodist(paper)
            P_bin = [P_cost_distribution(paper, c) for c in 0:m]
            emit("paper:empP", N_paper, P_emp)
            emit("paper:binP", N_paper, P_bin)

            println("done: $fam_name n=$n inst=$inst m=$m ceil=$(round(ceiling_ar; digits=4))")
            flush(stdout)
        end
    end

    suffix = tid == 0 ? (SMOKE ? "_smoke" : "_all") : "_task$tid"
    outpath = joinpath(@__DIR__, "results$suffix.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,c_opt,ceil_p1,variant," *
                    "mse_rel,ew_avg,ew_star_raw,ew_star_nrm," *
                    "ar_raw,regret_raw,dgamma_raw,dbeta_raw," *
                    "ar_nrm,regret_nrm,dgamma_nrm,dbeta_nrm," *
                    "pred_raw,normw_raw,fit_obj,fit_params")
        foreach(r -> println(io, r), rows)
    end
    println("E014 task(s) complete → $outpath")
end

main()
