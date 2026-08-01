#=
E047 — Harden the E046 polynomial-time N estimator: hyperparameter
robustness, larger-n validation, timing at scale.

Question:
  A. Is E046's headline (poly N + estimated counts + normalized objective is
     best at p=1, trails samp10 by ~0.005 at p=3) robust to the estimator
     hyperparameters? Sweep (S, K) in {(3,50),(10,50),(10,200),(10,800),
     (30,200)} at n=12. Does more profile sampling (K=800) close the p=3
     gap, or (per the E042 beneficial-noise story) hurt at p=1?
  B. Does the E044-style size trend hold at n=16 (poly's p=1 advantage grew
     from n=12 to n=14)? exact / samp10 / poly(S=10,K=200), both objectives,
     both depths, statevector truth still computable.
  C. Timing at scale: n in {16,20,24} x 3 ER(0.5) instances, poly estimator
     wall-clock plus Wang-Landau step/stage counts, vs the enumeration
     sampler where feasible (n<=20). No statevector at n=24.

Design: 7 families x 5 instances per part (seeds SEED=20260611,
instance_seed(fi,n,inst)=SEED+10000 fi+100 n+inst, matching E041/E044/E046);
p=1 grid (40x40) and p=3 linear-ramp grid (8^4); regret vs the statevector
grid ceiling. The (S=10,K=200) poly variant uses rng seed+999, exactly
reproducing E046's estimator draws. estimator.jl is include()d from E046,
not copied.

Run:   JULIA_NUM_THREADS=auto julia --project research/experiments/047_estimator-hardening/run.jl  (~45-70 min)
Smoke: E47_SMOKE=1 JULIA_NUM_THREADS=auto julia --project research/experiments/047_estimator-hardening/run.jl  (~3 min)

Outputs: results.csv (parts A and B), timing.csv (part C).
=#

using JuliaQAOA
using Random: MersenneTwister, Xoshiro
using Base.Threads: @threads

include(joinpath(@__DIR__, "..", "046_polytime-N-estimator", "estimator.jl"))

const SMOKE = get(ENV, "E47_SMOKE", "0") == "1"
const SEED  = 20260611
const INST  = SMOKE ? 2 : 5
const GLEN  = SMOKE ? 10 : 40
const RLEN  = SMOKE ? 4 : 8
const WLKW  = SMOKE ? (; lnf_final=1e-2, stage_cap=200_000) : (; lnf_final=1e-3, stage_cap=2_000_000)

# Part A: hyperparameter sweep (S, K, rng offset). Offset 999 for (10,200)
# reproduces E046's poly draws bit-for-bit at matching (n, instance).
const SK_FULL = [(3, 50, 901), (10, 50, 902), (10, 200, 999), (10, 800, 904), (30, 200, 905)]
const A_N  = SMOKE ? 10 : 12
const A_SK = SMOKE ? SK_FULL[[1, 3]] : SK_FULL
# Part B: larger-n validation with E046's best setting.
const B_N  = SMOKE ? 12 : 16   # smoke size must be even (3-regular family)
const B_SK = [(10, 200, 999)]
# Part C: timing scaling run.
const C_NS   = SMOKE ? [12] : [16, 20, 24]
const C_INST = SMOKE ? 1 : 3
const C_ENUM_MAX = 20     # enumeration sampler infeasible-by-policy above this

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

const P1γ = collect(range(0.0, π;   length = GLEN))
const P1β = collect(range(0.0, π/2; length = GLEN))
const Rγ  = collect(range(0.05, 1.6; length = RLEN))
const Rβ  = collect(range(0.05, 0.8; length = RLEN))
const P1_sched = vec([([g], [b]) for g in P1γ, b in P1β])
const R_combos = vec([(g1,gf,b1,bf) for g1 in Rγ, gf in Rγ, b1 in Rβ, bf in Rβ])
const R_sched  = [linear_ramp(c..., 3) for c in R_combos]

function objectives(N, γmat, βmat, counts::AbstractVector, copt)
    Q = QAOA_proxy_multi(N, γmat, βmat)[end]
    K = size(Q, 2)
    raw = zeros(Float64, K); nrm = zeros(Float64, K)
    for k in 1:K, ci in axes(Q, 1)
        w = counts[ci] * abs2(Q[ci, k])
        raw[k] += w * (ci - 1); nrm[k] += w
    end
    return raw ./ copt, raw ./ nrm ./ copt
end

"""
Prepare one instance: exact N, samp10 N (exact counts), and one poly variant
(estimated counts) per (S, K) setting. Fail fast if any poly run claims an
unattainable cost class; merely missing tiny classes is recorded.
"""
function prepare(fi, n, inst, sk_list)
    fam = FAMILIES[fi][1]
    seed = instance_seed(fi, n, inst)
    edges = FAMILIES[fi][2](MersenneTwister(seed), n)
    m = length(edges)
    @assert m > 0 "empty edge list: $fam n=$n inst=$inst"
    costs = maxcut_costs(n, edges)
    copt = maximum(costs)
    @assert copt > 0
    counts = zeros(Float64, m + 1)
    for c in costs; counts[Int(c)+1] += 1.0; end
    attained = findall(>(0.0), counts)
    t_ex  = @elapsed Nex = get_homogeneous_distribution_from_costs_direct(costs, m, n)
    t_s10 = @elapsed s10 = sampled_homogeneous_distribution(costs, m, n;
        samples_per_class=10, rng=MersenneTwister(seed + 777))
    variants = [("exact", Nex, counts, t_ex, 0), ("samp10", s10, counts, t_s10, 0)]
    for (S, K, off) in sk_list
        t = @elapsed begin
            Np, cp, _ = polytime_homogeneous_distribution(n, edges;
                samples_per_class=S, K_subsets=K, rng=Xoshiro(seed + off), WLKW...)
        end
        visited = findall(ci -> cp[ci] > 0.0, 1:(m + 1))
        @assert issubset(visited, attained) "poly estimator visited an unattainable cost class ($fam n=$n inst=$inst S=$S K=$K)"
        missed = length(setdiff(attained, visited))
        push!(variants, ("poly_S$(S)_K$(K)", Np, cp, t, missed))
    end
    (; fam, n, inst, m, costs, copt, variants)
end

"Run one validation part (A or B): all families x instances at size n."
function run_part(part::String, n::Int, sk_list)
    jobs = [(fi, inst) for fi in eachindex(FAMILIES) for inst in 1:INST]
    prep = map(jobs) do (fi, inst)
        q = prepare(fi, n, inst, sk_list)
        ts = join(("$(v[1])=$(round(v[4], digits=2))s" for v in q.variants), " ")
        println("prep $part: $(q.fam) n=$n inst=$inst m=$(q.m) $ts")
        q
    end
    rows = Vector{Vector{String}}(undef, length(prep))
    @threads for k in eachindex(prep)
        q = prep[k]
        out = String[]
        for (p, scheds, γmat, βmat) in (
                (1, P1_sched,
                    reshape([s[1][1] for s in P1_sched], :, 1),
                    reshape([s[2][1] for s in P1_sched], :, 1)),
                (3, R_sched,
                    linear_ramp_matrix([c[1] for c in R_combos], [c[2] for c in R_combos],
                                       [c[3] for c in R_combos], [c[4] for c in R_combos], 3)...))
            F = [qaoa_expectation(q.costs, q.n, γs, βs) / q.copt for (γs, βs) in scheds]
            ceilv = maximum(F)
            @assert isfinite(ceilv) && ceilv > 0
            for (kind, N, counts, t_est, missed) in q.variants
                fraw, fnorm = objectives(N, γmat, βmat, counts, q.copt)
                push!(out, join((q.fam, q.n, q.inst, part, p, round(ceilv, digits=6),
                    kind, round(t_est, digits=3), missed,
                    round(ceilv - F[argmax(fraw)], digits=6),
                    round(ceilv - F[argmax(fnorm)], digits=6)), ","))
            end
        end
        rows[k] = out
        println("done $part: ", q.fam, " n=", n, " inst=", q.inst)
    end
    return rows
end

"""
Mirror of estimator.jl's `wang_landau`, step for step and rand-call for
rand-call (the same seed yields the identical lng), adding step and stage
counters so the scaling run can report how much MCMC work the estimator
does. Kept here rather than in E046 so the frozen estimator stays untouched.
"""
function wang_landau_counted(n::Int, edges::Vector{Tuple{Int,Int}}; rng,
                             flat::Float64=0.8, lnf_final::Float64=1e-3,
                             check_every::Int=10_000 * n, stage_cap::Int=2_000_000)
    m = length(edges)
    adj = adjacency(n, edges)
    lng = fill(-Inf, m + 1)
    hist = zeros(Int, m + 1)
    s = rand(rng, 0:(1 << n) - 1)
    c = cut_cost(s, edges)
    lng[c + 1] = 0.0
    lnf = 1.0
    total_steps = 0
    stages = 0
    while lnf > lnf_final
        stages += 1
        steps = 0
        while steps < stage_cap
            for _ in 1:check_every
                v = rand(rng, 0:n-1)
                c2 = c + flip_delta(s, v, adj)
                lg2 = lng[c2 + 1]
                if isinf(lg2) || log(rand(rng)) < lng[c + 1] - lg2
                    s ⊻= 1 << v
                    c = c2
                    if isinf(lng[c + 1]); lng[c + 1] = 0.0; end
                end
                lng[c + 1] += lnf
                hist[c + 1] += 1
            end
            steps += check_every
            vis = findall(isfinite, lng)
            hv = view(hist, vis)
            if minimum(hv) >= flat * (sum(hv) / length(vis))
                break
            end
        end
        total_steps += steps
        fill!(hist, 0)
        lnf /= 2
    end
    return lng, total_steps, stages
end

function run_timing(io)
    for n in C_NS, inst in 1:C_INST
        seed = instance_seed(1, n, inst)
        edges = erdos_renyi_edges(n, 0.5; rng=MersenneTwister(seed))
        m = length(edges)
        t_enum = NaN
        attained = Int[]
        if n <= C_ENUM_MAX
            t_enum = @elapsed begin
                costs = maxcut_costs(n, edges)
                sampled_homogeneous_distribution(costs, m, n;
                    samples_per_class=10, rng=MersenneTwister(seed + 777))
            end
            attained = findall(c -> any(==(c - 1), costs), 1:(m + 1))
        end
        t_poly = @elapsed begin
            _, cpoly, _ = polytime_homogeneous_distribution(n, edges;
                samples_per_class=10, K_subsets=200, rng=Xoshiro(seed + 999), WLKW...)
        end
        # Counted rerun of the Wang-Landau stage only, identical rng stream.
        t_wl = @elapsed begin
            lng, wl_steps, wl_stages = wang_landau_counted(n, edges;
                rng=Xoshiro(seed + 999), WLKW...)
        end
        missed = ""
        if n <= C_ENUM_MAX
            visited = findall(ci -> cpoly[ci] > 0.0, 1:(m + 1))
            @assert issubset(visited, attained) "poly estimator visited an unattainable cost class (timing n=$n inst=$inst)"
            missed = string(length(setdiff(attained, visited)))
        end
        println(io, join(("ER(0.5)", n, inst, m,
            isnan(t_enum) ? "" : round(t_enum, digits=3),
            round(t_poly, digits=3), round(t_wl, digits=3),
            wl_steps, wl_stages, missed), ","))
        println("timing: n=$n inst=$inst m=$m t_enum=$(isnan(t_enum) ? "n/a" : round(t_enum, digits=2))s ",
                "t_poly=$(round(t_poly, digits=2))s t_wl=$(round(t_wl, digits=2))s ",
                "wl_steps=$wl_steps wl_stages=$wl_stages")
    end
end

function main()
    rows_a = run_part("A", A_N, A_SK)
    rows_b = run_part("B", B_N, B_SK)
    respath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(respath, "w") do io
        println(io, "family,n,inst,part,p,ceiling,kind,t_est_s,classes_missed,regret_raw,regret_norm")
        for rs in (rows_a, rows_b), r in rs, line in r; println(io, line); end
    end
    println("wrote $respath")

    timpath = joinpath(@__DIR__, SMOKE ? "timing_smoke.csv" : "timing.csv")
    open(timpath, "w") do io
        println(io, "family,n,inst,m,t_enum_s,t_poly_s,t_wl_s,wl_steps,wl_stages,classes_missed")
        run_timing(io)
    end
    println("wrote $timpath")
end

main()
