#=
E036 — From state error to regret: a two-point leakage certificate.

The sharpest referee objection (journal_readiness gap 3): the paper's error
calculus bounds STATE error while parameter quality is an argmax-transfer
(parameter-space) matter, so the theory "bounds the wrong error." This
experiment verifies the missing link, a rigorous two-point bound:

With F(θ) the true normalized objective, F̂(θ) the exact-compression proxy
objective computed from the NORMALIZED compressed state χ = φ/‖φ‖, θ* the
true argmax and θ̂ the proxy argmax over the same candidate set,

  regret = F(θ*) − F(θ̂) ≤ ε(θ*) + ε(θ̂)        [F̂(θ*) ≤ F̂(θ̂) kills the middle term]

  ε(θ) = |F(θ) − F̂(θ)| ≤ ‖ψ−χ‖ · (σ_ψ + sqrt(σ_χ² + δ²)) / c_opt,

where σ_ψ, σ_χ are the cost standard deviations in ψ and χ, δ = ⟨C⟩_ψ−⟨C⟩_χ,
and ‖ψ−χ‖ ≤ ‖ψ−φ‖ + (1−‖φ‖) ≤ 2 Σ_ℓ λ_ℓ(θ) is controlled by the leakage
calculus (Theorem 2). Proof: insert a = ⟨C⟩_ψ into
⟨ψ|C|ψ⟩−⟨χ|C|χ⟩ = ⟨ψ−χ|(C−a)|ψ⟩ + ⟨χ|(C−a)|ψ−χ⟩ and Cauchy-Schwarz.

So the calculus DOES bound regret, needing leakage at only two points (both
of which sit in the low-leakage corner in practice). The experiment measures
the bound's tightness at p=1 and p=3 across the standard 140 instances, and
asserts the inequality on every instance (a machine check of the proof).
Also reported: the crude version (σ's replaced by m), and how often the
normalized proxy argmax differs from the unnormalized one used elsewhere.

Run: JULIA_NUM_THREADS=auto julia --project research/experiments/036_regret-certificate/run.jl  (~10 min)
Smoke: E20_SMOKE=1 ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads

const SMOKE = get(ENV, "E20_SMOKE", "0") == "1"
const SEED  = 20260611
const NS    = SMOKE ? [10] : [12, 14]
const INST  = SMOKE ? 2 : 10
const GLEN  = SMOKE ? 10 : 40
const RLEN  = SMOKE ? 4 : 8

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

"Normalized proxy objective for every schedule column of a proxy-multi run."
function proxy_objectives(N, γmat, βmat, counts)
    Qs = QAOA_proxy_multi(N, γmat, βmat)
    Q = Qs[end]                                     # (num_costs, K)
    K = size(Q, 2)
    raw = zeros(Float64, K); nrm = zeros(Float64, K)
    for k in 1:K, ci in axes(Q, 1)
        w = counts[ci] * abs2(Q[ci, k])
        raw[k] += w * (ci - 1)
        nrm[k] += w
    end
    return raw ./ nrm, raw ./ (2.0^0)               # (normalized ⟨C⟩, unused)
end

"Cost mean and std of a statevector."
function cost_stats_state(ψ, costs)
    m1 = 0.0; m2 = 0.0
    for i in eachindex(ψ)
        w = abs2(ψ[i]); c = costs[i]
        m1 += w * c; m2 += w * c^2
    end
    return m1, sqrt(max(m2 - m1^2, 0.0))
end

"Cost mean and std of normalized compressed amplitudes Q with class sizes."
function cost_stats_Q(Q, counts)
    nrm = 0.0; m1 = 0.0; m2 = 0.0
    for ci in eachindex(Q)
        w = counts[ci] * abs2(Q[ci]); c = ci - 1
        nrm += w; m1 += w * c; m2 += w * c^2
    end
    m1 /= nrm; m2 /= nrm
    return m1, sqrt(max(m2 - m1^2, 0.0))
end

"ε(θ): certified upper bound on |F − F̂|(θ), plus the crude version."
function epsilon_at(costs, n, γs, βs, copt, m)
    traj = compressed_qaoa_trajectory(costs, n, γs, βs)
    dist_chi = traj.distance[end] + (1.0 - traj.compressed_norm[end])
    ψ = qaoa_statevector(costs, n, γs, βs)
    μψ, σψ = cost_stats_state(ψ, costs)
    μχ, σχ = cost_stats_Q(traj.Qs[end], traj.counts)
    δ = μψ - μχ
    eps = dist_chi * (σψ + sqrt(σχ^2 + δ^2)) / copt
    eps_crude = dist_chi * 2m / copt
    eps_leak = (2 * sum(traj.leakage)) * (σψ + sqrt(σχ^2 + δ^2)) / copt
    return eps, eps_crude, eps_leak
end

function instance_rows(fam, n, inst, m, costs, copt, N, counts)
    rows = String[]
    for (p, scheds, γmat, βmat) in (
            (1, P1_sched,
                reshape([s[1][1] for s in P1_sched], :, 1),
                reshape([s[2][1] for s in P1_sched], :, 1)),
            (3, R_sched,
                linear_ramp_matrix([c[1] for c in R_combos], [c[2] for c in R_combos],
                                   [c[3] for c in R_combos], [c[4] for c in R_combos], 3)...))
        Ltrue = [qaoa_expectation(costs, n, γs, βs) / copt for (γs, βs) in scheds]
        istar = argmax(Ltrue)
        Fhat_norm, _ = proxy_objectives(N, γmat, βmat, counts)
        ihat = argmax(Fhat_norm)
        # unnormalized argmax (the convention of E002/E014) for comparison
        Qs = QAOA_proxy_multi(N, γmat, βmat)
        vraw = vec(expectation(Qs[end], counts ./ (1 << n), n))
        ihat_raw = argmax(vraw)

        regret = Ltrue[istar] - Ltrue[ihat]
        e_star, ec_star, el_star = epsilon_at(costs, n, scheds[istar]..., copt, m)
        e_hat,  ec_hat,  el_hat  = epsilon_at(costs, n, scheds[ihat]...,  copt, m)
        bound = e_star + e_hat
        crude = ec_star + ec_hat
        leakb = el_star + el_hat
        # the theorem, machine-checked
        regret <= bound + 1e-9 || error("certificate violated: $fam n=$n inst=$inst p=$p " *
                                        "regret=$regret bound=$bound")
        push!(rows, join((fam, n, inst, m, p,
            round(regret, digits=6), round(bound, digits=6),
            round(bound / max(regret, 1e-12), digits=2),
            round(crude, digits=6), round(leakb, digits=6),
            Int(ihat != ihat_raw)), ","))
    end
    return rows
end

function main()
    jobs = [(fi, fam, n, inst) for (fi, (fam, _)) in enumerate(FAMILIES)
            for n in NS for inst in 1:INST]
    prep = map(jobs) do (fi, fam, n, inst)
        edges = FAMILIES[fi][2](MersenneTwister(instance_seed(fi, n, inst)), n)
        m = length(edges); costs = maxcut_costs(n, edges); copt = maximum(costs)
        N = get_homogeneous_distribution_from_costs_direct(costs, m, n)
        counts = zeros(Int, m + 1); for c in costs; counts[Int(c)+1] += 1; end
        (; fam, n, inst, m, costs, copt, N, counts)
    end
    println("prep done ($(length(prep)) instances)")
    out = Vector{Vector{String}}(undef, length(prep))
    @threads for k in eachindex(prep)
        q = prep[k]
        out[k] = instance_rows(q.fam, q.n, q.inst, q.m, q.costs, q.copt, q.N, q.counts)
        println("done: ", q.fam, " n=", q.n, " inst=", q.inst)
    end
    path = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,m,p,regret,bound,ratio,crude_bound,leakage_bound,argmax_moved")
        for rs in out, r in rs; println(io, r); end
    end
    println("wrote $path; certificate held on every instance and depth")
end

main()
