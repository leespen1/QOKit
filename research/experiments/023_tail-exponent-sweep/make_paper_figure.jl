#=
Paper figure for the tail-exponent boundary (E023): regret vs Pareto α at
p=1, both families. The binned proxy is nearly α-independent; best-of-3
rescaled transfer sweeps through it, crossing at the infinite-variance
point α = 2.
Run from repo root: julia --project research/experiments/023_tail-exponent-sweep/make_paper_figure.jl
=#
using Statistics: mean
using CairoMakie

const DIR = @__DIR__
reg = Dict{Tuple{String, Float64, String, Int}, Vector{Float64}}()   # (fam, α, method, inst)
for f in filter(f -> occursin(r"^results_task\d+\.csv$", f), readdir(DIR))
    for (i, line) in enumerate(eachline(joinpath(DIR, f)))
        i == 1 && continue
        p = split(line, ",")
        p[7] == "ceiling" && continue
        parse(Int, p[8]) == 1 || continue
        α = parse(Float64, replace(String(p[6]), "a" => ""))
        key = (String(p[1]), α, String(p[7]), parse(Int, p[3]))
        push!(get!(reg, key, Float64[]), parse(Float64, p[11]))
    end
end
αs = [1.2, 1.5, 2.0, 3.0, 5.0]
fams = ["ER(0.5)", "3-regular"]

"per-instance best-of-3 transfer, then cell mean"
function best_transfer(fam, α)
    vals = Float64[]
    for inst in 1:10
        srcs = [reg[(fam, α, "transfer_$r", inst)][1] for r in 1:3
                if haskey(reg, (fam, α, "transfer_$r", inst))]
        isempty(srcs) || push!(vals, minimum(srcs))
    end
    mean(vals)
end
cellmean(fam, α, m) = mean(vcat([reg[(fam, α, m, i)] for i in 1:10
                                 if haskey(reg, (fam, α, m, i))]...))

fig = Figure(size=(560, 360), fontsize=12)
ax = Axis(fig[1, 1]; xlabel="Pareto tail exponent α", ylabel="regret (p = 1)",
          xscale=log10, xticks=(αs, string.(αs)))
colors = CairoMakie.Makie.wong_colors()
for (k, fam) in enumerate(fams)
    proxy = [cellmean(fam, α, "binned_sampled") for α in αs]
    trans = [best_transfer(fam, α) for α in αs]
    scatterlines!(ax, αs, proxy; color=colors[k], marker=:circle,
                  label="$fam — binned proxy")
    scatterlines!(ax, αs, trans; color=colors[k], marker=:utriangle,
                  linestyle=:dash, label="$fam — transfer (best of 3)")
end
vlines!(ax, [2.0]; color=:gray55, linestyle=:dot)
text!(ax, 2.05, 0.105; text="infinite variance ← α = 2 → finite", color=:gray45,
      fontsize=11)
axislegend(ax; position=:rt, framevisible=false)
save(joinpath(DIR, "alpha_crossover.png"), fig; px_per_unit=3)
println("wrote alpha_crossover.png")
