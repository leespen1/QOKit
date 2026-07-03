#=
Paper figure for E3.2 (depth scaling). Two panels from results_task{1,2,3}.csv
(n = 16, 18, 20):

  Left:  family-pooled mean accumulated leakage Σλ vs layer at n=20, one curve
         per ramp, with the small-angle prediction Σ_ℓ β_ℓ γ_ℓ² (Theorem 3 /
         E2.1 law, scaled to the small-ramp endpoint) — within-run shape
         follows the schedule profile at working ramps; norm-loss saturation
         takes over at large ones.
  Right: family-pooled mean overlap vs layer for small and moderate ramps at
         n = 16, 18, 20 — fidelity at depth stays usable in the working regime.

Writes depth_scaling.png next to this script (copy to the paper repo's
Figures/generated/ when updating the manuscript).

Run from the repo root:
  julia --project research/experiments/011_depth-scaling/make_figure.jl
=#

using Statistics: mean
using CairoMakie

const DIR = @__DIR__
const TASKS = Dict(16 => "results_task1.csv", 18 => "results_task2.csv",
                   20 => "results_task3.csv")
const P = 30

# rows[(n, ramp, layer)] -> vector of (cum_leakage, overlap) over family×instance
function load()
    acc = Dict{Tuple{Int,String,Int},Vector{Tuple{Float64,Float64}}}()
    for (n, file) in TASKS
        for line in readlines(joinpath(DIR, file))[2:end]
            f = split(line, ',')
            ramp, layer = String(f[6]), parse(Int, f[7])
            cum, ovl = parse(Float64, f[9]), parse(Float64, f[11])
            push!(get!(acc, (n, ramp, layer), Tuple{Float64,Float64}[]), (cum, ovl))
        end
    end
    return acc
end

function main()
    acc = load()
    counts = unique(length(v) for v in values(acc))
    @assert counts == [70] "expected 7 families × 10 instances per cell, got $counts"

    fig = Figure(size = (1150, 430), fontsize = 16)

    ax1 = Axis(fig[1, 1]; xlabel = "layer ℓ", ylabel = "accumulated leakage Σλ",
               title = "n = 20, by ramp size")
    rampcolors = [("small", :navy), ("moderate", :dodgerblue),
                  ("large", :darkorange), ("extreme", :firebrick)]
    for (ramp, color) in rampcolors
        ys = [mean(first.(acc[(20, ramp, ℓ)])) for ℓ in 1:P]
        lines!(ax1, 1:P, ys; color, linewidth = 2.5, label = ramp)
    end
    # small-angle prediction for the small ramp: λ_ℓ ∝ β_ℓ γ_ℓ², scaled to the
    # observed endpoint (γ: 0.1→0.4, β: 0.4→0.1 — run.jl RAMPS)
    γℓ = [0.1 + (0.4 - 0.1) * ℓ / P for ℓ in 1:P]
    βℓ = [0.4 + (0.1 - 0.4) * ℓ / P for ℓ in 1:P]
    pred = cumsum(βℓ .* γℓ .^ 2)
    pred .*= mean(first.(acc[(20, "small", P)])) / pred[end]
    lines!(ax1, 1:P, pred; color = :gray, linestyle = :dash,
           linewidth = 1.8, label = "βγ² profile (Thm 3)")
    axislegend(ax1; position = :lt)

    ax2 = Axis(fig[1, 2]; xlabel = "layer ℓ", ylabel = "overlap |⟨ψ|φ⟩|²",
               title = "small (solid) and moderate (dashed) ramps")
    ncolors = [(16, :seagreen), (18, :steelblue), (20, :purple)]
    for (n, color) in ncolors
        for (ramp, style) in [("small", :solid), ("moderate", :dash)]
            ys = [mean(last.(acc[(n, ramp, ℓ)])) for ℓ in 1:P]
            lines!(ax2, 1:P, ys; color, linestyle = style, linewidth = 2.5,
                   label = ramp == "small" ? "n = $n" : nothing)
        end
    end
    ylims!(ax2, 0, 1.02)
    axislegend(ax2; position = :lb)

    out = joinpath(DIR, "depth_scaling.png")
    save(out, fig)
    println("wrote $out")
end

main()
