#=
E015 follow-up: derivation check for the cost-conditioning correction
Var(E[T|c]) (see research/v2_conditioning_derivation.md).

Key exact identity: with S = Σ_{(j,k)∈E} s_j s_k = m − 2c and
X_e = s_j s_k for edge e = (j,k),

    S² = m + Σ_{e≠f share a vertex} X_e X_f + Σ_{e≠f disjoint} X_e X_f
       = m + (T − 2m) + R₄            (shared-vertex pairs ARE the codegree form)
    ⇒  T = S² + m − R₄,   R₄ = Σ_{ordered disjoint edge pairs} X_e X_f.

So E[T|c] has an EXACT quadratic part (m−2c)²; the crude linear regression
(= 36τ²/m) misses it. Exact moments (uniform spins, monomial cancellation):

    Var(S) = m                E[S³] = 6τ            (τ = #triangles)
    Var(S²) = 2m(m−1) + 24q   (q = #4-cycles)
    Cov(T,S) = 6τ             Cov(T,S²) = Var(T) = 4P₁ + 16q,  P₁ = Σ_v C(deg_v,2)

Quadratic projection of E[T|c] onto span{1, S, S²} (any graph):

    Vq = (36τ²·D − 72τ²·W + m·W²) / (m·D − 36τ²),   W = Var(T), D = 2m(m−1)+24q.

This script verifies every moment identity as an EXACT integer identity, then
measures (by full 2^n enumeration): Var(E[T|c]), the linear part 36τ²/m, the
quadratic projection (empirical vs closed form Vq), cubic/quartic projections,
and the remainder. Run:
  julia --project research/experiments/015_v2-density-law/conditioning_derivation_check.jl
=#

using JuliaQAOA
using Random: MersenneTwister
using Statistics: mean
using Printf: @sprintf

const NS = [12, 13, 14]
const NINST = 5
const SEED = 20260705
const FAMILIES = [
    ("ER(0.5)",  (rng, n) -> erdos_renyi_edges(n, 0.5; rng)),
    ("ER(0.25)", (rng, n) -> erdos_renyi_edges(n, 0.25; rng)),
    ("BA(k=4)",  (rng, n) -> barabasi_albert_edges(n, 4; rng)),
]

# subgraph counts entering the closed forms
function graph_counts(edges, n)
    A = zeros(Int, n, n)
    for (a, b) in edges
        A[a+1, b+1] = 1
        A[b+1, a+1] = 1
    end
    m = length(edges)
    deg = vec(sum(A; dims = 2))
    C = A * A                                     # C[j,k] = codegree (j≠k)
    τ = sum(C .* A) ÷ 6                           # tr(A³)/6
    q = sum(binomial(C[j, k], 2) for j in 1:n for k in (j+1):n) ÷ 2   # #4-cycles
    P1 = sum(binomial(d, 2) for d in deg)         # #edge pairs sharing a vertex
    W = 2 * sum(j == k ? 0 : C[j, k]^2 for j in 1:n, k in 1:n)        # Var(T)
    # Φ = Σ_{j<k} A_jk · P3_jk, P3_jk = #simple 3-paths j→k (enters Cov(T,S³))
    A3 = C * A
    Φ = sum(C[j, k] * (A3[j, k] - A[j, k] * (deg[j] + deg[k] - 1))
            for j in 1:n for k in (j+1):n)
    return (; A, m, τ, q, P1, W, Φ)
end

# full enumeration: per-S-class counts/sums of T, plus exact integer moments
function enumerate_classes(A, n, m)
    N = 1 << n
    cnt = zeros(Int, 2m + 1)                      # index S + m + 1
    sumT = zeros(Int, 2m + 1)
    ΣT = 0; ΣT2 = 0; ΣS = 0; ΣS2 = 0; ΣS3 = 0; ΣS4 = 0; ΣTS = 0; ΣTS2 = 0; ΣTS3 = 0
    s = Vector{Int}(undef, n)
    for x in 0:(N-1)
        for i in 1:n
            s[i] = ifelse(((x >> (i - 1)) & 1) == 1, -1, 1)
        end
        T = 0; twoS = 0
        for i in 1:n
            ai = 0
            @inbounds for j in 1:n
                ai += A[i, j] * s[j]
            end
            T += ai * ai
            twoS += ai * s[i]
        end
        S = twoS ÷ 2
        k = S + m + 1
        cnt[k] += 1; sumT[k] += T
        ΣT += T; ΣT2 += T * T; ΣS += S; ΣS2 += S * S; ΣS3 += S^3; ΣS4 += S^4
        ΣTS += T * S; ΣTS2 += T * S * S; ΣTS3 += T * S^3
    end
    return (; cnt, sumT, ΣT, ΣT2, ΣS, ΣS2, ΣS3, ΣS4, ΣTS, ΣTS2, ΣTS3, N)
end

# variance captured by weighted LS fit of class means on polynomials of S up to degree k
function proj_var(Svals, w, y, deg)
    z = Svals ./ max(maximum(abs, Svals), 1)      # scale for conditioning
    X = [z[i]^j for i in eachindex(z), j in 0:deg]
    G = X' * (w .* X)
    coef = G \ (X' * (w .* y))
    fit = X * coef
    μ = sum(w .* fit)
    return sum(w .* (fit .- μ) .^ 2)
end

function analyze(edges, n)
    g = graph_counts(edges, n)
    (; m, τ, q, P1, W, Φ) = g
    e = enumerate_classes(g.A, n, m)

    # ---- exact integer moment identities (fail-fast) ----
    @assert e.ΣT == 2m * e.N                       "E[T] = 2m failed"
    @assert e.ΣT2 - 4m^2 * e.N == W * e.N          "Var(T) = 4P₁+16q failed"
    @assert W == 4P1 + 16q                         "codegree ↔ P₁,q count failed"
    @assert e.ΣS == 0 && e.ΣS2 == m * e.N          "Var(S) = m failed"
    @assert e.ΣS3 == 6τ * e.N                      "E[S³] = 6τ failed"
    @assert e.ΣS4 == (3m^2 - 2m + 24q) * e.N       "E[S⁴] closed form failed"
    @assert e.ΣTS == 6τ * e.N                      "Cov(T,S) = 6τ failed"
    @assert e.ΣTS2 - 2m * m * e.N == W * e.N       "Cov(T,S²) = Var(T) failed"
    # Cov(T,S³) = 6τ(3m−2) + 12Φ  ⇒  E[TS³] = 2m·6τ + 6τ(3m−2) + 12Φ
    @assert e.ΣTS3 == (30m * τ - 12τ + 12Φ) * e.N  "Cov(T,S³) closed form failed"

    # ---- measured Var(E[T|c]) and polynomial projections ----
    ks = findall(>(0), e.cnt)
    w = e.cnt[ks] ./ e.N
    Svals = Float64.(ks .- (m + 1))
    tbar = e.sumT[ks] ./ e.cnt[ks]
    VB = sum(w .* (tbar .- 2m) .^ 2)               # exact Var(E[T|c])

    lin = 36τ^2 / m                                # linear projection, closed form
    D = 2m * (m - 1) + 24q
    Vq = (36τ^2 * D - 72τ^2 * W + m * W^2) / (m * D - 36τ^2)  # quadratic, closed form
    p2 = proj_var(Svals, w, tbar, 2)
    @assert isapprox(p2, Vq; rtol = 1e-8) "empirical deg-2 projection ≠ closed form Vq"
    p1 = proj_var(Svals, w, tbar, 1)
    @assert isapprox(p1, lin; rtol = 1e-8) "empirical deg-1 projection ≠ 36τ²/m"
    p3 = proj_var(Svals, w, tbar, 3)
    p4 = proj_var(Svals, w, tbar, 4)
    return (; m, τ, q, P1, W, VB, lin, Vq, p3, p4)
end

function main()
    κ(p) = 18 * (1 / 2 + 2p - p^2) / (1 / 4 + 3p^2 / 2 - p^3)   # ER n→∞ constant
    println("ER asymptotic constants κ(p) = Vq·m/τ² :  κ(0.5) = $(round(κ(0.5), digits=2)),",
            "  κ(0.25) = $(round(κ(0.25), digits=2))   (crude linear = 36)")
    println()
    println(rpad("family", 10), rpad("n", 4), rpad("m", 5), rpad("τ", 5),
            rpad("VarE[T|c]", 11), rpad("lin=36τ²/m", 12), rpad("Vq(closed)", 12),
            rpad("deg4 proj", 11), rpad("Vq/VB", 7), rpad("deg4/VB", 8), "VB·m/τ²")
    agg = Dict{String,Vector{NTuple{4,Float64}}}()
    for (fi, (fam, gen)) in enumerate(FAMILIES), n in NS, inst in 1:NINST
        edges = gen(MersenneTwister(SEED + 10_000fi + 100n + inst), n)
        r = analyze(edges, n)
        ratio = r.τ > 0 ? r.VB * r.m / r.τ^2 : NaN
        println(rpad(fam, 10), rpad(n, 4), rpad(r.m, 5), rpad(r.τ, 5),
                rpad(@sprintf("%.1f", r.VB), 11), rpad(@sprintf("%.1f", r.lin), 12),
                rpad(@sprintf("%.1f", r.Vq), 12), rpad(@sprintf("%.1f", r.p4), 11),
                rpad(@sprintf("%.3f", r.Vq / r.VB), 7),
                rpad(@sprintf("%.3f", r.p4 / r.VB), 8),
                r.τ > 0 ? @sprintf("%.1f", ratio) : "--")
        push!(get!(agg, fam, NTuple{4,Float64}[]),
              (r.Vq / r.VB, r.p4 / r.VB, ratio, r.Vq * r.m / max(r.τ, 1)^2))
    end
    println("\nFamily means:")
    for (fam, _) in FAMILIES
        v = agg[fam]
        println("  ", rpad(fam, 10),
                " Vq/VB = ", @sprintf("%.3f", mean(first.(v))),
                "   deg4/VB = ", @sprintf("%.3f", mean(getindex.(v, 2))),
                "   VB·m/τ² = ", @sprintf("%.1f", mean(filter(!isnan, getindex.(v, 3)))),
                "   Vq·m/τ² = ", @sprintf("%.1f", mean(getindex.(v, 4))))
    end
    println("\nAll exact integer moment identities held on every instance.")
end

main()
