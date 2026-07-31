#=
Paper figure for the fitted-shape/model-error story (§5.4): at MATCHED
entrywise error, regret is set by the error's direction relative to the
dynamics, not its size. One panel: mean normalized-objective regret vs the
perturbation scale ε for the four direction profiles of E017 Part A, pooled
over all 150 instances (bars = standard error), with the exact-N baseline.

Output: fitted_shape_directions.png (copy to the paper repo's
Figures/generated/).

Run from repo root: julia --project research/experiments/017_error-directions/make_paper_figure.jl
=#

using Statistics: mean, std
using CairoMakie

const DIR = @__DIR__
const PROFILES = ["lowd", "random", "highd", "nulled"]
const LABELS = Dict("lowd" => "aligned (low-d)", "random" => "random",
                    "highd" => "high-d (d≈n/2)", "nulled" => "nulled (invisible)")
const EPS = [0.01, 0.05, 0.2, 0.5]

regrets = Dict((s, e) => Float64[] for s in PROFILES, e in EPS)
exactn = Float64[]
for f in filter(f -> occursin(r"^results_task\d+\.csv$", f), readdir(DIR))
    for (i, line) in enumerate(eachline(joinpath(DIR, f)))
        i == 1 && continue
        p = split(line, ",")
        v = String(p[8])
        reg = parse(Float64, p[18])   # regret_nrm
        if v == "exactN"
            push!(exactn, reg)
        elseif startswith(v, "pert:")
            _, shape, eps = split(v, ":")
            push!(regrets[(String(shape), parse(Float64, eps))], reg)
        end
    end
end
@assert length(exactn) == 150 "expected 150 instances, got $(length(exactn))"

fig = Figure(size=(520, 340), fontsize=12)
ax = Axis(fig[1, 1];
          xlabel="entrywise perturbation scale  ε  (‖ΔN‖ = ε‖N‖)",
          ylabel="regret (normalized objective)",
          xscale=log10, xticks=(EPS, string.(EPS)))
colors = CairoMakie.Makie.wong_colors()
for (k, s) in enumerate(PROFILES)
    μ = [mean(regrets[(s, e)]) for e in EPS]
    se = [std(regrets[(s, e)]) / sqrt(length(regrets[(s, e)])) for e in EPS]
    errorbars!(ax, EPS, μ, se; color=colors[k], whiskerwidth=6)
    scatterlines!(ax, EPS, μ; color=colors[k], label=LABELS[s], marker=:circle)
end
hlines!(ax, [mean(exactn)]; color=:gray40, linestyle=:dash)
text!(ax, 0.011, mean(exactn) + 0.004; text="exact-N baseline",
      color=:gray40, fontsize=10)
axislegend(ax; position=:lt, framevisible=false, rowgap=0)
save(joinpath(DIR, "fitted_shape_directions.png"), fig; px_per_unit=3)
println("wrote $(joinpath(DIR, "fitted_shape_directions.png"))")
save(joinpath(DIR, "fitted_shape_directions.pdf"), fig)
println("wrote $(joinpath(DIR, "fitted_shape_directions.pdf"))")
