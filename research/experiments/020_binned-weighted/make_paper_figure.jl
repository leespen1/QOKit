#=
Paper figure for the weighted-MaxCut binning extension (E020). Two panels:
(a) one-layer leakage vs K at three angle pairs (ER(0.5), n=16 — the
    hardest cell): the structural floor plus the O(1/K) binning variance;
(b) binned-proxy regret vs K for all four cells at p=1 (solid) and p=3
    (dashed), with the unweighted exact-compression reference band.

Output: binned_weighted.png (copy to the paper repo's Figures/generated/).
Run from repo root: julia --project research/experiments/020_binned-weighted/make_paper_figure.jl
=#

using Statistics: mean
using CairoMakie

const DIR = @__DIR__
leak = Dict{Tuple{String, Int, String, Int}, Vector{Float64}}()
reg = Dict{Tuple{String, Int, Int, Int}, Vector{Float64}}()
for f in filter(f -> occursin(r"^results_task\d+\.csv$", f), readdir(DIR))
    for (i, line) in enumerate(eachline(joinpath(DIR, f)))
        i == 1 && continue
        p = split(line, ",")
        fam, n, kind, K = String(p[1]), parse(Int, p[2]), String(p[6]), parse(Int, p[7])
        if startswith(kind, "leak")
            push!(get!(leak, (fam, n, kind, K), Float64[]), parse(Float64, p[9]))
        elseif kind == "binned_proxy"
            push!(get!(reg, (fam, n, parse(Int, p[8]), K), Float64[]), parse(Float64, p[11]))
        end
    end
end
const KS = [4, 8, 16, 32, 64, 128]

fig = Figure(size=(980, 340), fontsize=12)
colors = CairoMakie.Makie.wong_colors()

ax1 = Axis(fig[1, 1]; xlabel="number of bins K", ylabel="one-layer leakage λ",
           xscale=log2, yscale=log10, xticks=(KS, string.(KS)),
           title="(a) leakage vs K   —   G(16, 0.5), U[0,1] weights")
angles = [("leak_g0.2_b0.2", "γ=0.2, β=0.2"), ("leak_g0.5_b0.3", "γ=0.5, β=0.3"),
          ("leak_g1.0_b0.4", "γ=1.0, β=0.4")]
for (k, (kind, lab)) in enumerate(angles)
    μ = [mean(leak[("ER(0.5)", 16, kind, K)]) for K in KS]
    scatterlines!(ax1, KS, μ; color=colors[k], label=lab, marker=:circle)
end
ref = [0.09 / sqrt(K / 4) for K in KS]
lines!(ax1, KS, ref; color=:gray55, linestyle=:dot)
text!(ax1, 8.5, 0.028; text="∝ 1/√K", color=:gray45, fontsize=12)
axislegend(ax1; position=:lb, framevisible=false)

ax2 = Axis(fig[1, 2]; xlabel="number of bins K", ylabel="regret",
           xscale=log2, xticks=(KS, string.(KS)),
           title="(b) binned-proxy regret vs K")
cells = [("ER(0.5)", 14), ("ER(0.5)", 16), ("3-regular", 14), ("3-regular", 16)]
for (k, (fam, n)) in enumerate(cells)
    μ1 = [mean(reg[(fam, n, 1, K)]) for K in KS]
    μ3 = [mean(reg[(fam, n, 3, K)]) for K in KS]
    scatterlines!(ax2, KS, μ1; color=colors[k], label="$fam n=$n", marker=:circle)
    scatterlines!(ax2, KS, μ3; color=colors[k], linestyle=:dash, marker=:utriangle)
end
hspan!(ax2, 0.028, 0.055; color=(:gray, 0.15))
text!(ax2, 4.2, 0.058; text="unweighted exact-compression range (p=1)",
      color=:gray40, fontsize=10)
axislegend(ax2; position=:rt, framevisible=false)
text!(ax2, 24, 0.135; text="dashed: p=3 ramps\nsolid: p=1", color=:gray40,
      fontsize=11)

save(joinpath(DIR, "binned_weighted.png"), fig; px_per_unit=3)
println("wrote binned_weighted.png")
