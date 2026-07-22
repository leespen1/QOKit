#=
E021 — Does the section-4 theory survive weighted MaxCut?

The paper's Limits scope everything to unweighted MaxCut. The theory chain
(neighbor-sum lemma -> variance reduction -> codegree form -> quadratic
conditioning bound -> cubic leakage law) should generalize verbatim to
INTEGER edge weights with the replacements
  m   -> W        (total weight)
  m2  -> sum w^2, q4 -> sum w^4
  A_jk -> A^w_jk = sum_i w_ij w_ik      (weighted codegree)
  tau -> tau_w = sum_triangles w_ab w_bc w_ca
  c4  -> c4_w  = sum_4cycles  w_ab w_bc w_cd w_da
i.e.
  (L1w)  sum_i c(x xor e_i) = (n-4) c(x) + 2W
  (L2w)  s2 = (n-8) c^2 + 4 W c + T,  V2 = within-class Var(T)
  (L3w)  T = sum_{j,k} A^w_jk s_j s_k;  E[T] = 2 m2;  Var(T) = 2 sum_{j!=k} (A^w_jk)^2
  (M)    Var(S) = m2;  E[S^3] = 6 tau_w;  Var(S^2) = 2 m2^2 - 2 q4 + 24 c4_w;
         Cov(T,S) = 6 tau_w;  Cov(T,S^2) = Var(T)
  (P)    Var(E[T|c]) >= v' G^{-1} v with v = (6 tau_w, Var T),
         G = [[m2, 6 tau_w], [6 tau_w, Var(S^2)]]
  (C)    lambda_1 ~ (|beta| gamma^2 / 8) sqrt(V2) at small angles.

Integer weights keep the cost classes exponentially large; CONTINUOUS
weights collapse them to near-singletons and dissolve the compression
itself (a structural boundary, checked here by measuring the class-size
profile for uniform(0,1) weights).

All identities are asserted (1e-9); the bound capture and cubic-law ratio
are measured across seven families with iid weights in {1,2}, n = 12, 14.

Run: JULIA_NUM_THREADS=auto julia --project research/experiments/021_weighted-maxcut/run.jl  (~5 min)
=#

using JuliaQAOA
using Random: MersenneTwister, rand
using Base.Threads: @threads
using LinearAlgebra: tr

const SEED = 20260611
const NS   = [12, 14]
const INST = 10
const WCHOICES = [1.0, 2.0]

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

relerr(a, b) = abs(a - b) / max(abs(a), abs(b), 1e-12)

"Weighted graph quantities from the weight matrix Wm (n x n, symmetric)."
function graph_quantities_w(Wm)
    n = size(Wm, 1)
    W  = sum(Wm) / 2
    m2 = sum(abs2, Wm) / 2
    q4 = sum(x -> x^4, Wm) / 2
    B2 = Wm * Wm                       # (B2)_{jk} = A^w_{jk} incl. diagonal
    tau_w = tr(B2 * Wm) / 6
    # weighted 4-cycles: tr(Wm^4) counts closed 4-walks with weight products:
    # tr(W^4) = 8 c4_w + 2 sum_{j!=k} A^w_jk * ...; safer: use the identity
    # tr(W^4) = sum_{jk} (B2_{jk})^2, and remove non-cycle walks:
    # closed 4-walks = pairs of 2-walks: sum_{jk} B2_{jk}^2; degenerate walks
    # are (i) j=k contributions where the two 2-walks pass any middle vertex:
    # handled below exactly as in the unweighted case with weight powers.
    trW4 = sum(abs2, B2)
    offB2sq = trW4 - sum(abs2, [B2[i, i] for i in 1:n])   # sum_{j!=k} (A^w_jk)^2
    # c4_w from the standard formula generalized: 8 c4_w = tr(W^4) - sum_j (A^w_jj)^2
    #   - sum_{j!=k} [ 2 * Wm[j,k]^2 * ... ]  -- for weighted graphs the clean
    # route is direct enumeration over vertex quadruples (n <= 14: fine).
    c4_w = 0.0
    for a in 1:n, b in 1:n, c in 1:n, d in 1:n
        if a < b && a < c && a < d && b != c && c != d && b != d
            # count each unordered 4-cycle once: fix a as minimum, orient b-c-d
            # cycle a-b-c-d-a; each undirected cycle counted twice (two dirs)
            c4_w += Wm[a, b] * Wm[b, c] * Wm[c, d] * Wm[d, a] / 2
        end
    end
    varT_formula = 2.0 * offB2sq
    return (; W, m2, q4, tau_w, c4_w, varT_formula)
end

"Exact enumeration: cost classes, moments of (S, T), Var(E[T|c]), V2."
function exact_w(n, edges, w)
    nbrw = [Tuple{Int,Float64}[] for _ in 1:n]
    for (e, (i, j)) in enumerate(edges)
        push!(nbrw[i + 1], (j + 1, w[e])); push!(nbrw[j + 1], (i + 1, w[e]))
    end
    W = sum(w)
    # integer classes: costs are multiples of 1 (weights in {1,2})
    key(c) = round(Int, 2c)            # robust integer key (half-steps safe)
    acc = Dict{Int,Vector{Float64}}()  # key -> [count, sumT, sumT2]
    ES = 0.0; ES2 = 0.0; ES3 = 0.0; ES4 = 0.0
    ET = 0.0; ET2 = 0.0; ETS = 0.0; ETS2 = 0.0
    σ = Vector{Float64}(undef, n)
    neighbor_lemma_checked = false
    for y in 0:(2^n - 1)
        for i in 1:n
            σ[i] = 1.0 - 2.0 * ((y >> (i - 1)) & 1)
        end
        S = 0.0
        for (e, (i, j)) in enumerate(edges)
            S += w[e] * σ[i + 1] * σ[j + 1]
        end
        c = (W - S) / 2
        T = 0.0
        for i in 1:n
            s = 0.0
            for (j, wij) in nbrw[i]
                s += wij * σ[j]
            end
            T += s * s
        end
        if !neighbor_lemma_checked   # spot-check L1w on the first bitstring
            nb = 0.0
            for i in 1:n
                Si = S
                for (j, wij) in nbrw[i]
                    Si -= 2 * wij * σ[i] * σ[j]
                end
                nb += (W - Si) / 2
            end
            relerr(nb, (n - 4) * c + 2W) < 1e-9 || error("L1w fails")
            neighbor_lemma_checked = true
        end
        v = get!(acc, key(c), zeros(3))
        v[1] += 1; v[2] += T; v[3] += T^2
        ES += S; ES2 += S^2; ES3 += S^3; ES4 += S^4
        ET += T; ET2 += T^2; ETS += T * S; ETS2 += T * S^2
    end
    N = 2.0^n
    ES /= N; ES2 /= N; ES3 /= N; ES4 /= N
    ET /= N; ET2 /= N; ETS /= N; ETS2 /= N
    varT = ET2 - ET^2
    varE = 0.0; V2 = 0.0; nclasses = 0; maxclass = 0
    for (_, v) in acc
        nclasses += 1; maxclass = max(maxclass, Int(v[1]))
        μ = v[2] / v[1]
        varE += v[1] / N * (μ - ET)^2
        V2 += v[1] / N * (v[3] / v[1] - μ^2)
    end
    return (; W, ES, ES2, ES3, ES4, ET, ETS, ETS2, varT, varE, V2, nclasses, maxclass)
end

function main()
    jobs = [(fi, fam, n, inst) for (fi, (fam, _)) in enumerate(FAMILIES)
            for n in NS for inst in 1:INST]
    rows = Vector{String}(undef, length(jobs))
    @threads for k in eachindex(jobs)
        fi, fam, n, inst = jobs[k]
        rng = MersenneTwister(instance_seed(fi, n, inst))
        edges = FAMILIES[fi][2](rng, n)
        w = [WCHOICES[rand(rng, 1:2)] for _ in edges]
        Wm = zeros(n, n)
        for (e, (i, j)) in enumerate(edges)
            Wm[i + 1, j + 1] = w[e]; Wm[j + 1, i + 1] = w[e]
        end
        G = graph_quantities_w(Wm)
        M = exact_w(n, edges, w)

        # identities
        relerr(M.ET, 2 * G.m2) < 1e-9 || error("E[T]=2m2 fails: $fam $n $inst")
        relerr(M.varT, G.varT_formula) < 1e-9 || error("weighted codegree Var(T) fails: $fam $n $inst")
        relerr(M.ES2, G.m2) < 1e-9 || error("Var(S)=m2 fails: $fam $n $inst")
        relerr(M.ES3, 6 * G.tau_w) < 1e-9 || error("E[S^3]=6tau_w fails: $fam $n $inst")
        relerr(M.ES4 - M.ES2^2, 2 * G.m2^2 - 2 * G.q4 + 24 * G.c4_w) < 1e-9 ||
            error("Var(S^2) weighted 4-cycle identity fails: $fam $n $inst")
        relerr(M.ETS - M.ET * M.ES, 6 * G.tau_w) < 1e-9 || error("Cov(T,S)=6tau_w fails: $fam $n $inst")
        relerr(M.ETS2 - M.ET * M.ES2, M.varT) < 1e-9 || error("Cov(T,S^2)=Var(T) fails: $fam $n $inst")
        # law of total variance closes: Var(T) = V2 + Var(E[T|c])
        relerr(M.varT, M.V2 + M.varE) < 1e-9 || error("total variance fails: $fam $n $inst")

        # quadratic bound capture
        g11 = G.m2; g12 = 6 * G.tau_w
        g22 = 2 * G.m2^2 - 2 * G.q4 + 24 * G.c4_w
        v1 = 6 * G.tau_w; v2 = M.varT
        qb = (g22 * v1^2 - 2 * g12 * v1 * v2 + g11 * v2^2) / (g11 * g22 - g12^2)
        qb <= M.varE + 1e-6 || error("bound violated: $fam $n $inst qb=$qb varE=$(M.varE)")
        capture = qb / M.varE

        # cubic law at small angles: lambda_1 vs (|b| g^2 / 8) sqrt(V2)
        γ0, β0 = 0.05, 0.05
        costs = Float64[]
        # build integer-keyed cost vector for the trajectory (costs are
        # integers or half-integers times 1; weights in {1,2} give integers)
        cvec = zeros(Float64, 1 << n)
        for y in 0:(2^n - 1)
            c = 0.0
            for (e, (i, j)) in enumerate(edges)
                bi = (y >> i) & 1; bj = (y >> j) & 1
                c += w[e] * (bi != bj ? 1.0 : 0.0)
            end
            cvec[y + 1] = c
        end
        traj = compressed_qaoa_trajectory(cvec, n, [γ0], [β0])
        λ1 = traj.leakage[1]
        pred = abs(β0) * γ0^2 / 8 * sqrt(M.V2)
        cubic_ratio = λ1 / pred

        rows[k] = join((fam, n, inst, round(G.W, digits=1), M.nclasses, M.maxclass,
            round(M.varT, digits=4), round(M.varE, digits=4), round(capture, digits=4),
            round(M.V2, digits=4), round(cubic_ratio, digits=4)), ",")
        println("done: $fam n=$n inst=$inst capture=$(round(capture, digits=3)) cubic=$(round(cubic_ratio, digits=3))")
    end
    path = joinpath(@__DIR__, "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,W,nclasses,maxclass,varT,varE_exact,capture_quad,V2,cubic_ratio")
        for r in rows; println(io, r); end
    end

    # continuous-weight degeneracy check (one instance per family, n=12)
    println("\ncontinuous-weight class collapse (uniform(0,1) weights, n=12):")
    open(joinpath(@__DIR__, "continuous_weights.csv"), "w") do io
        println(io, "family,nclasses_over_2n,maxclass")
        for (fi, (fam, gen)) in enumerate(FAMILIES)
            rng = MersenneTwister(instance_seed(fi, 12, 1))
            edges = gen(rng, 12)
            wu = [rand(rng) for _ in edges]
            seen = Dict{Float64,Int}()
            for y in 0:(2^12 - 1)
                c = 0.0
                for (e, (i, j)) in enumerate(edges)
                    bi = (y >> i) & 1; bj = (y >> j) & 1
                    c += wu[e] * (bi != bj ? 1.0 : 0.0)
                end
                seen[round(c; digits=10)] = get(seen, round(c; digits=10), 0) + 1
            end
            frac = length(seen) / 2^12
            println(io, "$fam,$(round(frac, digits=4)),$(maximum(values(seen)))")
            println("  $fam: distinct classes = $(round(100 * frac, digits=1))% of 2^n, largest class = $(maximum(values(seen)))")
        end
    end
    println("wrote results.csv and continuous_weights.csv; all weighted identities held")
end

main()
