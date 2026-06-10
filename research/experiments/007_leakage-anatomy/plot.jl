#=
Figure for experiment 007: (left) one-layer leakage λ_uniform over the (γ, β)
plane for ER(0.5) at n=14 (instance mean), with the empirical-proxy argmax and
the analytical proxy's spurious peak (exp 004) marked; (right) small-angle
scaling: λ_uniform vs γ at fixed small β for all families, log-log, showing
the slope-2 law (λ ~ β·γ²) and the density ordering.

Run from the repo root: julia --project research/experiments/007_leakage-anatomy/plot.jl
=#

using CairoMakie
using Statistics: mean

dir = @__DIR__
lines_csv = readlines(joinpath(dir, "results.csv"))
hdr = split(lines_csv[1], ',')
rows = [split(l, ',') for l in lines_csv[2:end]]
col(name) = findfirst(==(name), hdr)
gf(r, name) = parse(Float64, r[col(name)])

γvals = sort(unique(gf(r, "gamma") for r in rows))
βvals = sort(unique(gf(r, "beta") for r in rows))
fams = unique(r[1] for r in rows)

fig = Figure(size=(1150, 460))

# Left: ER(0.5) n=14 heatmap of mean λ_uniform
er = [r for r in rows if r[1] == "ER(0.5)" && r[2] == "14"]
λmap = [mean(gf(r, "lambda_uniform") for r in er
             if gf(r, "gamma") == γ && gf(r, "beta") == β) for γ in γvals, β in βvals]
ax1 = Axis(fig[1, 1]; xlabel="γ", ylabel="β",
           title="One-layer leakage λ(γ, β) — ER(0.5), n = 14 (mean of 10 instances)")
hm = heatmap!(ax1, γvals, βvals, λmap)
Colorbar(fig[1, 2], hm)
scatter!(ax1, [0.403], [0.282]; color=:cyan, marker=:star5, markersize=18,
         label="exact-compression argmax")
scatter!(ax1, [0.644], [1.168]; color=:red, marker=:xcross, markersize=16,
         label="analytical proxy's spurious peak")
axislegend(ax1; position=:rt, labelsize=10, backgroundcolor=(:white, 0.7))

# Right: small-angle scaling, λ vs γ at fixed small β, all families (n=14)
β0 = βvals[3]
ax2 = Axis(fig[1, 3]; xlabel="γ", ylabel="λ_uniform at β ≈ $(round(β0, digits=2))",
           xscale=log10, yscale=log10,
           title="Small-angle scaling: slope ≈ 2 (λ ~ β·γ²)")
for (i, f) in enumerate(sort(fams))
    sel = [r for r in rows if r[1] == f && r[2] == "14" && gf(r, "beta") == β0]
    λs = [mean(gf(r, "lambda_uniform") for r in sel if gf(r, "gamma") == γ) for γ in γvals[1:10]]
    lines!(ax2, γvals[1:10], λs; label=f, color=Makie.wong_colors()[mod1(i, 7)])
end
# slope-2 guide
g = γvals[1:6]
lines!(ax2, g, 0.8 .* (g ./ g[1]) .^ 2 .* 1e-3; color=:black, linestyle=:dash, label="slope 2")
axislegend(ax2; position=:rb, labelsize=9)

save(joinpath(dir, "leakage_anatomy.png"), fig; px_per_unit=2)
println("saved → ", joinpath(dir, "leakage_anatomy.png"))
