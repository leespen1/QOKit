#=
E045 — A selection-aware sharpening of the two-point regret certificate.

Question: E036's certificate regret <= eps(theta*) + eps(theta-hat) is
~10x loose; E038's exact anatomy regret = [e(theta*) - e(theta-hat)] - margin
shows a nonnegative margin = F-hat(theta-hat) - F-hat(theta*) that the proxy
itself reports. Subtracting it is free and rigorous:

  regret = [e(theta*) - e(theta-hat)] - margin            (exact identity)
        <= |e(theta*)| + |e(theta-hat)| - margin          (drop signs only)
        <= eps(theta*) + eps(theta-hat) - margin          (eps >= |e| pointwise)

where eps(theta) is E036's Cauchy-Schwarz bound
  eps(theta) = (||psi-phi|| + 1 - ||phi||) (sigma_psi + sqrt(sigma_chi^2 + delta^2)) / c_opt.
The margin subtraction uses no inequality, so the sharpened bound is valid
whenever the paper's is, and never worse. Two additional variants measured:
(a) the practitioner-uniform bound eps(theta-hat) + max_theta [eps(theta)
    - margin(theta)], which removes the oracle theta* entirely (rigorous:
    instantiate the per-theta inequality F(theta) - F(theta-hat) <=
    eps(theta) + eps(theta-hat) - margin(theta) at theta = theta*, then
    take the sup over the candidate set);
(b) the non-certified one-sided diagnostic eps(theta-hat) - margin, valid
    exactly when e(theta*) <= 0 (overprediction at theta*, which E038 found
    at 277/280 points but which no proxy-side quantity certifies).
Also measured: the cosine alignment of the two error vectors
psi - chi at theta* and theta-hat, to gauge headroom for a correlated-error
(two-point Cauchy-Schwarz) refinement.

This run asserts eps >= |e| at EVERY schedule of every landscape (a
pointwise machine check of the eps lemma), the exact identity, margin >= 0,
and validity of the paper, sharpened, and uniform bounds on every row.

Run: JULIA_NUM_THREADS=auto julia --project research/experiments/045_selection-aware-certificate/run.jl (~15 min)
Smoke: E45_SMOKE=1 ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads
using Statistics: median, mean

const SMOKE = get(ENV, "E45_SMOKE", "0") == "1"
const SEED  = 20260611
const NS    = SMOKE ? [10] : [12, 14]
const INST  = SMOKE ? 2 : 5
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

"True F, F-hat (normalized compressed), and eps(theta) for one schedule."
function point_quantities(costs, n, γs, βs, copt)
    ψ = qaoa_statevector(costs, n, γs, βs)
    μψ = 0.0; m2ψ = 0.0
    for i in eachindex(ψ)
        w = abs2(ψ[i]); c = costs[i]
        μψ += w * c; m2ψ += w * c^2
    end
    σψ = sqrt(max(m2ψ - μψ^2, 0.0))

    traj = compressed_qaoa_trajectory(costs, n, γs, βs)
    Q, counts = traj.Qs[end], traj.counts
    nrm = 0.0; μχ = 0.0; m2χ = 0.0
    for ci in eachindex(Q)
        w = counts[ci] * abs2(Q[ci]); c = ci - 1
        nrm += w; μχ += w * c; m2χ += w * c^2
    end
    μχ /= nrm; m2χ /= nrm
    σχ = sqrt(max(m2χ - μχ^2, 0.0))
    δ = μψ - μχ
    dist_chi = traj.distance[end] + (1.0 - traj.compressed_norm[end])
    eps = dist_chi * (σψ + sqrt(σχ^2 + δ^2)) / copt

    F    = μψ / copt
    Fhat = μχ / copt
    return F, Fhat, eps
end

"Cosine alignment of the error vectors psi - chi at two schedules."
function error_alignment(costs, n, sched_a, sched_b)
    errs = map((sched_a, sched_b)) do (γs, βs)
        ψ = qaoa_statevector(costs, n, γs, βs)
        traj = compressed_qaoa_trajectory(costs, n, γs, βs)
        Q = traj.Qs[end] ./ traj.compressed_norm[end]   # chi in compressed rep
        err = similar(ψ)
        for i in eachindex(ψ)
            err[i] = ψ[i] - Q[Int(costs[i]) + 1]
        end
        err
    end
    ip = sum(conj(a) * b for (a, b) in zip(errs[1], errs[2]))
    na = sqrt(sum(abs2, errs[1])); nb = sqrt(sum(abs2, errs[2]))
    return real(ip) / max(na * nb, 1e-300)
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
        K = length(scheds)
        F = zeros(K); Fhat = zeros(K); eps = zeros(K)
        for k in 1:K
            F[k], Fhat[k], eps[k] = point_quantities(costs, n, scheds[k]..., copt)
            # the eps lemma, machine-checked pointwise (e has arbitrary sign)
            abs(F[k] - Fhat[k]) <= eps[k] + 1e-10 ||
                error("eps < |e| at $fam n=$n inst=$inst p=$p k=$k")
        end
        # consistency with the proxy-multi landscape used by E036/E038
        Qend = QAOA_proxy_multi(N, γmat, βmat)[end]
        for k in 1:K
            raw = 0.0; nrm = 0.0
            for ci in axes(Qend, 1)
                w = counts[ci] * abs2(Qend[ci, k])
                raw += w * (ci - 1); nrm += w
            end
            abs(raw / nrm / copt - Fhat[k]) < 1e-8 ||
                error("proxy-multi vs trajectory F-hat mismatch: $fam n=$n inst=$inst p=$p k=$k")
        end

        istar = argmax(F); ihat = argmax(Fhat)
        regret = F[istar] - F[ihat]
        estar  = F[istar] - Fhat[istar]
        ehat   = F[ihat]  - Fhat[ihat]
        margin = Fhat[ihat] - Fhat[istar]
        margin >= -1e-12 || error("negative margin: $fam n=$n inst=$inst p=$p")
        abs(regret - ((estar - ehat) - margin)) < 1e-9 ||
            error("identity fails: $fam n=$n inst=$inst p=$p")

        bound_paper = eps[istar] + eps[ihat]
        bound_sharp = bound_paper - margin
        margins = Fhat[ihat] .- Fhat                     # margin(theta) >= 0 for all theta
        minimum(margins) >= -1e-12 || error("negative margin(theta): $fam n=$n inst=$inst p=$p")
        bound_unif = eps[ihat] + maximum(eps .- margins)

        # validity, machine-checked on every row
        regret <= bound_sharp + 1e-9 ||
            error("sharpened bound violated: $fam n=$n inst=$inst p=$p regret=$regret bound=$bound_sharp")
        bound_sharp <= bound_paper + 1e-12 || error("sharpened exceeds paper bound")
        bound_unif >= bound_sharp - 1e-9 ||
            error("uniform bound below two-point sharpened: $fam n=$n inst=$inst p=$p")

        # non-certified one-sided diagnostic (valid iff e(theta*) <= 0)
        onesided = eps[ihat] - margin
        onesided_sign_ok = Int(estar <= 0)
        onesided_holds = Int(regret <= onesided + 1e-9)

        errcos = error_alignment(costs, n, scheds[istar], scheds[ihat])

        push!(rows, join((fam, n, inst, m, p,
            round(regret, digits=6), round(bound_paper, digits=6),
            round(margin, digits=6), round(bound_sharp, digits=6),
            round(bound_unif, digits=6),
            round(bound_paper / max(regret, 1e-12), digits=2),
            round(bound_sharp / max(regret, 1e-12), digits=2),
            round(estar, digits=6), round(ehat, digits=6),
            onesided_sign_ok, round(onesided, digits=6), onesided_holds,
            round(errcos, digits=4)), ","))
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
    header = "family,n,inst,m,p,regret,bound_paper,margin,bound_sharp,bound_unif," *
             "ratio_paper,ratio_sharp,e_star,e_hat,onesided_sign_ok,onesided_bound," *
             "onesided_holds,err_cos"
    open(path, "w") do io
        println(io, header)
        for rs in out, r in rs; println(io, r); end
    end
    println("wrote $path; all validity assertions held")

    # summary
    cols = split(header, ",")
    data = [split(r, ",") for rs in out for r in rs]
    col(name) = [parse(Float64, d[findfirst(==(name), cols)]) for d in data]
    for p in (1, 3)
        sel = [parse(Int, d[findfirst(==("p"), cols)]) == p for d in data]
        reg = col("regret")[sel]
        pos = reg .>= 1e-4                        # ratio meaningful only off the regret floor
        rp = col("ratio_paper")[sel]; rs_ = col("ratio_sharp")[sel]
        bp = col("bound_paper")[sel]; bs = col("bound_sharp")[sel]
        bu = col("bound_unif")[sel];  mg = col("margin")[sel]
        oh = col("onesided_holds")[sel]; os = col("onesided_sign_ok")[sel]
        ec = col("err_cos")[sel]
        twice = count(bs .* 2 .<= bp) / length(bs)
        println("p=$p ($(count(pos))/$(length(reg)) rows with regret>=1e-4): ",
                "median ratio paper=$(round(median(rp[pos]), digits=2)) ",
                "sharp=$(round(median(rs_[pos]), digits=2)) ",
                "max sharp=$(round(maximum(rs_[pos]), digits=1)); ",
                "mean bound paper=$(round(mean(bp), digits=3)) ",
                "sharp=$(round(mean(bs), digits=3)) ",
                "unif=$(round(mean(bu), digits=3)); ",
                "mean margin=$(round(mean(mg), digits=3)); ",
                "frac >=2x tighter=$(round(twice, digits=3)); ",
                "onesided sign ok $(count(os .== 1))/$(length(os)), ",
                "holds $(count(oh .== 1))/$(length(oh)); ",
                "median err_cos=$(round(median(ec), digits=3))")
    end
end

main()
