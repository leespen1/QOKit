#=
E035 — Toward closing the density law: an exact lower bound on the
conditioning correction Var(E[T|c]).

E015 reduced the open V2 problem to one quantity: the conditioning
correction Var(E[T|c]) = Var(T) - V2 (law of total variance), measured to
track tau^2/m with a constant ~60 on dense families vs. the crude analytic
36 (linear-projection estimate). This experiment tests a sharper, still
fully rigorous estimate: project T onto span{S, S^2} where S = sum of edge
spin products = m - 2c. Since E[T|S] is the L2 projection of T onto ALL
functions of S, the polynomial projection variance is an exact LOWER bound:

  Var(E[T|c]) >= v' G^{-1} v,
  v = (Cov(T,S), Cov(T,S^2)) = (6 tau, Var_unc(T)),
  G = [[m, 6 tau], [6 tau, 2m^2 - 2m + 24 c4]],

with tau = triangles, c4 = 4-cycles, Var_unc(T) = 2 sum_{j != k} A_jk^2
(E015 Lemma). The three new moment identities used are themselves exact:
  Cov(T, S)   = 6 tau
  Cov(T, S^2) = Var_unc(T)          <- the surprise that makes this work
  E[S^3] = 6 tau,  Var(S^2) = 2m^2 - 2m + 24 c4.
All are asserted numerically per instance (1e-9 relative) before use.

Questions: (a) do the identities hold? (b) what fraction of the exact
Var(E[T|c]) does the quadratic bound capture, per family? (c) does the
bound explain why the measured constant exceeds the crude 36 (i.e., the
missing piece was the quadratic channel, present even in triangle-free
graphs)?

Run: julia --project research/experiments/035_conditioning-correction/run.jl  (~2 min)
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads
using LinearAlgebra: tr

const SEED = 20260611                        # same instances as E002/E004/E014/E034
const NS   = [12, 14]
const INST = 10

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

function graph_quantities(n, edges)
    B = zeros(Int, n, n)
    for (i, j) in edges          # edges are 0-indexed (QOKit convention)
        B[i + 1, j + 1] = 1; B[j + 1, i + 1] = 1
    end
    B2 = B * B
    m = length(edges)
    tau = tr(B2 * B) ÷ 6
    offB2 = sum(B2) - tr(B2)
    c4 = (tr(B2 * B2) - 2m - 2 * offB2) ÷ 8
    A2sum = sum(abs2, B2) - sum(abs2, [B2[i, i] for i in 1:n])  # sum_{j!=k} A_jk^2
    return m, tau, c4, A2sum
end

"Enumerate all 2^n strings: exact moments of (S, T) and Var(E[T|S])."
function exact_moments(n, edges)
    nbr = [Int[] for _ in 1:n]
    for (i, j) in edges          # 0-indexed edges
        push!(nbr[i + 1], j + 1); push!(nbr[j + 1], i + 1)
    end
    m = length(edges)
    # accumulators per S value (S in -m:m); index by cut count k = (m - S)/2 in 0:m
    cntS = zeros(Int, m + 1); sumT = zeros(Float64, m + 1)
    ES = 0.0; ES2 = 0.0; ES3 = 0.0; ES4 = 0.0
    ET = 0.0; ET2 = 0.0; ETS = 0.0; ETS2 = 0.0
    σ = Vector{Int}(undef, n)
    for y in 0:(2^n - 1)
        for i in 1:n
            σ[i] = 1 - 2 * ((y >> (i - 1)) & 1)
        end
        S = 0
        for (i, j) in edges
            S += σ[i + 1] * σ[j + 1]
        end
        T = 0
        for i in 1:n
            s = 0
            for j in nbr[i]
                s += σ[j]
            end
            T += s * s
        end
        k = (m - S) ÷ 2
        cntS[k + 1] += 1; sumT[k + 1] += T
        ES += S; ES2 += S^2; ES3 += S^3; ES4 += S^4
        ET += T; ET2 += T^2; ETS += T * S; ETS2 += T * S^2
    end
    N = 2.0^n
    ES /= N; ES2 /= N; ES3 /= N; ES4 /= N
    ET /= N; ET2 /= N; ETS /= N; ETS2 /= N
    varT = ET2 - ET^2
    varE = 0.0   # Var(E[T|S])
    for k in 0:m
        cntS[k + 1] == 0 && continue
        μ = sumT[k + 1] / cntS[k + 1]
        varE += cntS[k + 1] / N * (μ - ET)^2
    end
    return (; ES, ES2, ES3, ES4, ET, ETS, ETS2, varT, varE)
end

relerr(a, b) = abs(a - b) / max(abs(a), abs(b), 1e-12)

function main()
    jobs = [(fi, fam, n, inst) for (fi, (fam, _)) in enumerate(FAMILIES)
            for n in NS for inst in 1:INST]
    rows = Vector{String}(undef, length(jobs))
    @threads for k in eachindex(jobs)
        fi, fam, n, inst = jobs[k]
        edges = FAMILIES[fi][2](MersenneTwister(instance_seed(fi, n, inst)), n)
        m, tau, c4, A2sum = graph_quantities(n, edges)
        M = exact_moments(n, edges)

        # --- assert the exact identities (fail fast) ---
        varT_formula = 2.0 * A2sum
        relerr(M.varT, varT_formula) < 1e-9 ||
            error("Var(T) codegree identity fails: $fam n=$n inst=$inst")
        relerr(M.ET, 2.0 * m) < 1e-9 || error("E[T]=2m fails")
        relerr(M.ES2, Float64(m)) < 1e-9 || error("Var(S)=m fails")
        relerr(M.ES3, 6.0 * tau) < 1e-9 || error("E[S^3]=6tau fails: $fam n=$n inst=$inst")
        relerr(M.ES4 - M.ES2^2, 2.0 * m^2 - 2.0 * m + 24.0 * c4) < 1e-9 ||
            error("Var(S^2) 4-cycle identity fails: $fam n=$n inst=$inst")
        covTS = M.ETS - M.ET * M.ES
        relerr(covTS, 6.0 * tau) < 1e-9 || error("Cov(T,S)=6tau fails: $fam n=$n inst=$inst")
        covTS2 = M.ETS2 - M.ET * M.ES2
        relerr(covTS2, M.varT) < 1e-9 ||
            error("Cov(T,S^2)=Var(T) fails: $fam n=$n inst=$inst")

        # --- the bounds ---
        crude = m > 0 ? 36.0 * tau^2 / m : 0.0          # linear channel only
        g11 = Float64(m); g12 = 6.0 * tau
        g22 = 2.0 * m^2 - 2.0 * m + 24.0 * c4
        det = g11 * g22 - g12^2
        v1 = 6.0 * tau; v2 = M.varT
        quadbound = (g22 * v1^2 - 2 * g12 * v1 * v2 + g11 * v2^2) / det
        capture = quadbound / M.varE
        lincapture = (v1^2 / g11) / M.varE
        rho_cond_exact = 1.0 - M.varE / M.varT          # = V2 / Var(T)
        rho_cond_bound = 1.0 - quadbound / M.varT       # upper bound on V2/Var(T)

        rows[k] = join((fam, n, inst, m, tau, c4,
            round(M.varT, digits=4), round(M.varE, digits=4),
            round(v1^2 / g11, digits=4), round(crude, digits=4),
            round(quadbound, digits=4),
            round(lincapture, digits=4), round(capture, digits=4),
            round(rho_cond_exact, digits=4), round(rho_cond_bound, digits=4)), ",")
        println("done: $fam n=$n inst=$inst capture=$(round(capture, digits=3))")
    end
    path = joinpath(@__DIR__, "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,m,tau,c4,varT,varE_exact,linproj,crude36," *
                    "quadbound,capture_lin,capture_quad,rho_cond_exact,rho_cond_bound")
        for r in rows; println(io, r); end
    end
    println("wrote $path ($(length(rows)) instances); all identities held to 1e-9")
end

main()
