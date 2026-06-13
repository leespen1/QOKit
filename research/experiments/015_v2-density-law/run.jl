#=
E015 — Why does compression error scale with density? An exact decomposition of
the cubic-leakage variance functional V2, and the mechanism behind √V2 ∝ m.

Background. Corollary (cubic) gives the first-layer MaxCut leakage
  λ1 = (|β| γ² / 8) √V2 + h.o.,   V2 = 2^-n Σ_v M_v Var_{S_v}[s2],
with s2(y) = Σ_i c(y⊕e_i)². The paper measures √V2 ∝ m ("density, not
ER-ness") but leaves an analytic estimate of V2 open. This experiment proves two
exact reductions and uses them to explain the density law mechanistically.

Exact results verified here on every instance (fail-fast):
  L1 (reduction).  V2 = 2^-n Σ_v M_v Var_{S_v}[T],  T(y) = Σ_i δ_i(y)²,
      δ_i(y) = c(y⊕e_i) − c(y) = deg(i) − 2κ_i(y) = Σ_{j∼i} s_j  (s = ±1 spins),
      because s2(y) = (n−8)c(y)² + 4m c(y) + T(y) and the first two terms are
      class-constant (Lemma: Σ_i c(y⊕e_i) = (n−4)c(y) + 2m).
  L2 (codegree closed form).  Var_y(T) = 2 Σ_{j≠k} A_{jk}²  over UNIFORM y,
      where A_{jk} = #common neighbors of j,k. (Proof: T−2m = 2 Σ_{j<k}A_{jk}s_js_k.)
  C  (cubic law).  λ1 ≈ (|β|γ²/8)√V2 at small angles (measured ratio ≈ 1).

The mechanism (measured): V2 is the WITHIN-cost-class variance, i.e. a fraction
ρ_cond = V2 / Var_y(T) of the unconditional codegree variance. Var_y(T) grows
super-linearly with density (∝ density²·m² for ER), but ρ_cond SHRINKS with
density (more triangles ⇒ cost explains more of T), and the two effects cancel
to leave √V2 ≈ const·m across families — the density law.

Outputs results.csv with one row per instance:
  family,n,inst,seed,m,V2,V2_from_s2,varT_uncond,varT_codegree,SA2,rho_cond,
  lambda1,lambda1_pred,lambda_ratio

Run from the repo root (CPU only, ~1–2 min):
  julia --project research/experiments/015_v2-density-law/run.jl
Smoke test: E15_SMOKE=1 julia --project research/experiments/015_v2-density-law/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister
using Statistics: mean, var
using Printf: @sprintf

const SMOKE = get(ENV, "E15_SMOKE", "0") == "1"
const SEED  = 20260611                       # SAME instances as exp 002/004
const NS    = SMOKE ? [10, 12] : [12, 14, 16, 18]
const INST  = SMOKE ? 2 : 10
const SMALLγ, SMALLβ = 0.05, 0.05            # small-angle regime for the cubic law

const FAMILIES = [
    ("ER(0.5)",       (rng, n) -> erdos_renyi_edges(n, 0.5; rng)),
    ("ER(0.25)",      (rng, n) -> erdos_renyi_edges(n, 0.25; rng)),
    ("BA(k=2)",       (rng, n) -> barabasi_albert_edges(n, 2; rng)),
    ("BA(k=4)",       (rng, n) -> barabasi_albert_edges(n, 4; rng)),
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("WS(k=4;b=0.5)", (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]
instance_seed(fam_idx, n, inst) = SEED + 10_000 * fam_idx + 100 * n + inst

# O(2^n) within-class population variance of vals[y+1] grouped by integer cvec.
function within_class_var(vals::Vector{Float64}, cvec::Vector{Int}, n::Int)
    cmax = maximum(cvec)
    cnt = zeros(Int, cmax + 1); s = zeros(cmax + 1); s2 = zeros(cmax + 1)
    @inbounds for k in eachindex(vals)
        c = cvec[k] + 1; cnt[c] += 1; s[c] += vals[k]; s2[c] += vals[k]^2
    end
    tot = 0.0
    @inbounds for c in 1:cmax+1
        cnt[c] < 2 && continue
        μ = s[c] / cnt[c]
        tot += s2[c] - cnt[c] * μ^2          # M_c * Var_{S_c}
    end
    return tot / (1 << n)
end

function analyze(edges, n)
    costs = maxcut_costs(n, edges); m = length(edges)
    cy(y) = costs[y + 1]
    s2 = Vector{Float64}(undef, 1 << n)
    T  = Vector{Float64}(undef, 1 << n)
    cvec = Vector{Int}(undef, 1 << n)
    @inbounds for y in 0:(1 << n)-1
        cy0 = costs[y + 1]; acc2 = 0.0; accT = 0.0
        for i in 1:n
            ci = costs[(y ⊻ (1 << (i - 1))) + 1]
            acc2 += ci^2; accT += (ci - cy0)^2
        end
        s2[y + 1] = acc2; T[y + 1] = accT; cvec[y + 1] = Int(cy0)
    end
    V2       = within_class_var(T,  cvec, n)
    V2_s2    = within_class_var(s2, cvec, n)
    varT_unc = var(T; corrected = false)

    # codegree closed form: Var_y(T) = 2 Σ_{j≠k} A_jk²
    adj = [Set{Int}() for _ in 1:n]
    for e in edges; push!(adj[e[1] + 1], e[2] + 1); push!(adj[e[2] + 1], e[1] + 1); end
    SA2 = 0
    for j in 1:n, k in 1:n
        j == k && continue
        SA2 += length(intersect(adj[j], adj[k]))^2
    end
    varT_codeg = 2 * SA2

    # measured first-layer leakage at small angles
    tr = compressed_qaoa_trajectory(costs, n, [SMALLγ], [SMALLβ])
    λ1 = tr.leakage[1]
    λ1_pred = (abs(SMALLβ) * SMALLγ^2 / 8) * sqrt(V2)

    return (; m, V2, V2_s2, varT_unc, varT_codeg, SA2,
            rho = V2 / varT_unc, λ1, λ1_pred, ratio = λ1 / λ1_pred)
end

function main()
    rows = String[]
    push!(rows, "family,n,inst,seed,m,V2,V2_from_s2,varT_uncond,varT_codegree,SA2,rho_cond,lambda1,lambda1_pred,lambda_ratio")
    for (fi, (fam, gen)) in enumerate(FAMILIES), n in NS, inst in 1:INST
        seed = instance_seed(fi, n, inst)
        edges = gen(MersenneTwister(seed), n)
        r = analyze(edges, n)
        # --- fail-fast on the EXACT identities (machine precision) ---
        @assert isapprox(r.V2, r.V2_s2; rtol = 1e-10) "L1 (s2→T reduction) failed: $fam n=$n inst=$inst"
        @assert isapprox(r.varT_unc, r.varT_codeg; rtol = 1e-10) "L2 (codegree form) failed: $fam n=$n inst=$inst"
        @assert isapprox(r.λ1, r.λ1_pred; rtol = 0.03) "cubic law off by >3%: $fam n=$n inst=$inst ratio=$(r.ratio)"
        push!(rows, join((fam, n, inst, seed, r.m,
            @sprintf("%.6f", r.V2), @sprintf("%.6f", r.V2_s2),
            @sprintf("%.6f", r.varT_unc), @sprintf("%.6f", r.varT_codeg), r.SA2,
            @sprintf("%.6f", r.rho), @sprintf("%.8f", r.λ1),
            @sprintf("%.8f", r.λ1_pred), @sprintf("%.6f", r.ratio)), ","))
        println("done: ", rpad(fam, 14), " n=$n inst=$inst  m=$(r.m) V2=$(round(r.V2,digits=1)) ρ=$(round(r.rho,digits=3)) λ-ratio=$(round(r.ratio,digits=4))")
    end
    out = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(out, "w") do io; for r in rows; println(io, r); end; end
    println("\nwrote $out  ($(length(rows)-1) instances)")
    println("All exact identities (L1, L2) held to 1e-10; cubic law to 3% on every instance.")
end

main()
