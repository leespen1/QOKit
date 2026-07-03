#=
Two-panel figure for E015 (V₂ anatomy). Reads results.csv.

  (a) The cubic law: measured first-layer leakage λ₁ (at γ=β=0.05) vs the
      prediction (|β|γ²/8)√V₂, all 280 instances — they lie on y=x.
  (b) The density-law cancellation: per family (ordered by edge count), the
      UNCONDITIONAL √Var(T)/m (= √(2ΣA²)/m, codegree form) spreads ~2× across
      families, while the actual within-class √V₂/m is flattened to ~1.4× by the
      density-dependent conditioning survival ρ_cond — so √V₂ ≈ const·m.

Writes v2_anatomy.png next to this script (copy to the paper repo's
Figures/generated/ if used in the manuscript).

Run from the repo root:
  julia --project research/experiments/015_v2-density-law/make_figure.jl
=#

using DelimitedFiles: readdlm
using Statistics: mean
using CairoMakie

const DIR = @__DIR__

function main()
    M, hdr = readdlm(joinpath(DIR, "results.csv"), ','; header = true)
    cols = vec(hdr)
    col(name) = M[:, findfirst(==(name), cols)]
    fam = String.(col("family")); nn = Int.(col("n"))
    m = Float64.(col("m")); V2 = Float64.(col("V2")); varT = Float64.(col("varT_uncond"))
    λ1 = Float64.(col("lambda1")); λpred = Float64.(col("lambda1_pred"))

    fig = Figure(size = (1150, 460), fontsize = 15)

    # ---- (a) cubic law ----
    axa = Axis(fig[1, 1]; xlabel = "predicted  (|β|γ²/8)·√V₂",
        ylabel = "measured first-layer leakage λ₁",
        title = "(a) cubic law λ₁ = (|β|γ²/8)√V₂  (γ=β=0.05)")
    lo, hi = extrema(vcat(λ1, λpred))
    lines!(axa, [lo, hi], [lo, hi]; color = :gray60, linestyle = :dash, label = "y = x")
    scatter!(axa, λpred, λ1; color = (:navy, 0.5), markersize = 7)
    axislegend(axa; position = :lt)
    text!(axa, 0.04, 0.86; text = "280 instances; ratio ∈ [0.997, 1.002]",
          space = :relative, fontsize = 13)

    # ---- (b) flattening by family (n=16) ----
    N = 16
    FAMS = unique(fam)
    fam_m   = [mean(m[(fam .== f) .& (nn .== N)]) for f in FAMS]
    sqV_m   = [sqrt(mean(V2[(fam .== f) .& (nn .== N)]))   / fam_m[i] for (i, f) in enumerate(FAMS)]
    sqVarT_m= [sqrt(mean(varT[(fam .== f) .& (nn .== N)])) / fam_m[i] for (i, f) in enumerate(FAMS)]
    ord = sortperm(fam_m; rev = true)            # dense → sparse
    famlabel = Dict("ER(0.5)"=>"ER(0.5)","ER(0.25)"=>"ER(0.25)","BA(k=2)"=>"BA(2)",
        "BA(k=4)"=>"BA(4)","WS(k=4;b=0.1)"=>"WS(0.1)","WS(k=4;b=0.5)"=>"WS(0.5)",
        "3-regular"=>"3-reg")
    xs = 1:length(FAMS)
    axb = Axis(fig[1, 2]; ylabel = "√(leakage variance) per edge",
        title = "(b) conditioning flattens √leakage/m toward const (n=16)",
        xticks = (xs, [famlabel[FAMS[ord[i]]] for i in xs]),
        xticklabelrotation = π/5)
    # dumbbells
    for i in xs
        f = ord[i]
        lines!(axb, [i, i], [sqV_m[f], sqVarT_m[f]]; color = :gray70)
    end
    scatter!(axb, xs, [sqVarT_m[ord[i]] for i in xs]; color = :firebrick, markersize = 13,
             marker = :utriangle, label = "√Var(T)/m  (unconditional, 2ΣA²)")
    scatter!(axb, xs, [sqV_m[ord[i]] for i in xs]; color = :navy, markersize = 13,
             marker = :circle, label = "√V₂/m  (within-class; the density law)")
    axislegend(axb; position = :rt, labelsize = 12)
    ylims!(axb, 0, nothing)

    out = joinpath(DIR, "v2_anatomy.png")
    save(out, fig)
    println("wrote $out")
end

main()
