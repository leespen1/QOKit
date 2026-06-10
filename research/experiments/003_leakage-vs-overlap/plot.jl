#=
Figure for experiment 003: Theorem-2 bound tightness and per-layer behavior.
Left: accumulated leakage vs actual trajectory distance at p=20 (all 840 runs),
with the bound line y = x. Right: per-layer medians for the moderate ramp at
n=14, by family.

Run from the repo root: julia --project research/experiments/003_leakage-vs-overlap/plot.jl
=#

using CairoMakie
using Statistics: median

dir = @__DIR__
lines = readlines(joinpath(dir, "results.csv"))
hdr = split(lines[1], ',')
rows = [split(l, ',') for l in lines[2:end]]
col(name) = findfirst(==(name), hdr)
gf(r, name) = parse(Float64, r[col(name)])

fin = [r for r in rows if r[col("layer")] == "20"]
ramps = ["small", "moderate", "large", "extreme"]
ramp_colors = Dict(zip(ramps, Makie.wong_colors()[1:4]))

fig = Figure(size=(1100, 450))

ax1 = Axis(fig[1, 1]; xlabel="accumulated leakage Σλ_ℓ", ylabel="‖ψ₂₀ − φ₂₀‖",
           title="Theorem-2 bound at p = 20 (840 runs)")
for ramp in ramps
    sel = [r for r in fin if r[col("ramp")] == ramp]
    scatter!(ax1, [gf(r, "cum_leakage") for r in sel], [gf(r, "distance") for r in sel];
             color=ramp_colors[ramp], markersize=5, label=ramp)
end
maxx = maximum(gf(r, "cum_leakage") for r in fin)
lines!(ax1, [0, maxx], [0, maxx]; color=:black, linestyle=:dash, label="bound y = x")
axislegend(ax1; position=:rb)

ax2 = Axis(fig[1, 2]; xlabel="layer ℓ", ylabel="median ‖ψ_ℓ − φ_ℓ‖",
           title="Distance growth by family (moderate ramp, n = 14)")
fams = unique(r[1] for r in rows)
for (i, f) in enumerate(sort(fams))
    sel = [r for r in rows if r[1] == f && r[2] == "14" && r[col("ramp")] == "moderate"]
    med = [median(gf(r, "distance") for r in sel if r[col("layer")] == string(ℓ)) for ℓ in 1:20]
    lines!(ax2, 1:20, med; label=f, color=Makie.wong_colors()[mod1(i, 7)])
end
axislegend(ax2; position=:lt, labelsize=10)

save(joinpath(dir, "bound_tightness.png"), fig; px_per_unit=2)
println("saved → ", joinpath(dir, "bound_tightness.png"))
