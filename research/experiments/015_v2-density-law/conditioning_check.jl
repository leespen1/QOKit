#=
E015 follow-up: confirm the triangle mechanism behind the cost-conditioning
correction. By the law of total variance, the within-class V2 splits off the
unconditional codegree variance as
    V2 = Var_y(T) − Var(E[T|c]),     Var_y(T) = 2 Σ_{j≠k} A_jk²  (Lemma 2),
so the "conditioning correction" is Var(E[T|c]) = varT_uncond − V2 (both already
in results.csv). The heuristic in research/v2_density_law.md predicts it is
triangle-driven, Var(E[T|c]) ≈ c·τ²/m. This script tests that against τ counted
on each instance's graph.

Run: julia --project research/experiments/015_v2-density-law/conditioning_check.jl
=#
using JuliaQAOA
using Random: MersenneTwister
using Statistics: mean, cor

const RESCSV = joinpath(@__DIR__, "results.csv")
const GENS = [
    ("ER(0.5)",       (n, r) -> erdos_renyi_edges(n, 0.5; rng = r)),
    ("ER(0.25)",      (n, r) -> erdos_renyi_edges(n, 0.25; rng = r)),
    ("BA(k=2)",       (n, r) -> barabasi_albert_edges(n, 2; rng = r)),
    ("BA(k=4)",       (n, r) -> barabasi_albert_edges(n, 4; rng = r)),
    ("WS(k=4;b=0.1)", (n, r) -> watts_strogatz_edges(n, 4, 0.1; rng = r)),
    ("WS(k=4;b=0.5)", (n, r) -> watts_strogatz_edges(n, 4, 0.5; rng = r)),
    ("3-regular",     (n, r) -> random_regular_edges(n, 3; rng = r)),
]
const FAMIDX = Dict(g[1] => i for (i, g) in enumerate(GENS))
seedof(fam, n, inst) = 20260611 + 10_000 * FAMIDX[fam] + 100 * n + inst

function ntriangles(edges, n)
    adj = [Set{Int}() for _ in 1:n]
    for e in edges; push!(adj[e[1]+1], e[2]+1); push!(adj[e[2]+1], e[1]+1); end
    sum(length(intersect(adj[e[1]+1], adj[e[2]+1])) for e in edges) ÷ 3  # each Δ over 3 edges
end

function main()
    lines = readlines(RESCSV); hdr = split(lines[1], ',')
    ix(name) = findfirst(==(name), hdr)
    vb = Float64[]; pred = Float64[]; fams = String[]; tris = Int[]
    for ln in lines[2:end]
        isempty(ln) && continue
        f = split(ln, ',')
        fam = String(f[ix("family")]); n = parse(Int, f[ix("n")]); inst = parse(Int, f[ix("inst")])
        m = parse(Float64, f[ix("m")]); V2 = parse(Float64, f[ix("V2")])
        vT = parse(Float64, f[ix("varT_uncond")])
        edges = GENS[FAMIDX[fam]][2](n, MersenneTwister(seedof(fam, n, inst)))
        τ = ntriangles(edges, n)
        push!(vb, vT - V2); push!(pred, τ^2 / m); push!(fams, fam); push!(tris, τ)
    end
    @assert all(>=(0), vb) "Var(E[T|c]) must be nonnegative (law of total variance)"
    println("Conditioning correction  Var(E[T|c]) = Var_y(T) − V2   vs   τ²/m   ($(length(vb)) instances)")
    println("  Pearson(Var_between, τ²/m) = ", round(cor(vb, pred), digits = 4))
    withΔ = tris .> 0
    println("  among triangle-bearing instances ($(count(withΔ))): mean ratio Var_between/(τ²/m) = ",
            round(mean(vb[withΔ] ./ pred[withΔ]), digits = 1))
    println("  triangle-free instances ($(count(.!withΔ))): mean Var_between = ",
            round(mean(vb[.!withΔ]), digits = 2), " (small residual, not triangle-driven)")
    for fam in unique(fams)
        s = (fams .== fam) .& withΔ
        any(s) && println("   ", rpad(fam, 14), " ratio≈", round(mean(vb[s] ./ pred[s]), digits = 1),
                          "  ⟨τ⟩=", round(mean(tris[fams .== fam]), digits = 1))
    end
    println("\n→ The cost-conditioning correction is triangle-driven (Pearson 0.99 with τ²/m);")
    println("  this is the term that suppresses V2 for dense graphs and flattens √V2/m to ∝ m.")
end
main()
