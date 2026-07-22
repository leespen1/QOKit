#=
Paper figure for E2.4 (the fitted-shape paradox), for the model-error
subsection. Two panels over the 9 non-empirical models × 140 instances at p=1:

  Left:  regret vs the fit objective (slice-normalized entrywise MSE, log x) —
         no correlation; the fitted Triangle sits at lower MSE and higher
         regret than its unfitted default.
  Right: regret vs argmax displacement ‖(Δγ, Δβ)‖ from the empirical proxy's
         argmax — the quantity that actually predicts regret (ρ ≈ 0.7).

Reads results.csv and landscape.csv; writes fitted_shape_paradox.png next to
this script (copy to the paper repo's Figures/generated/ when updating the
manuscript).

Run from the repo root:
  julia --project research/experiments/012_fitted-shape-paradox/make_paper_figure.jl
=#

using CairoMakie

const DIR = @__DIR__

function read_csv(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    col = Dict(h => i for (i, h) in enumerate(header))
    return [split(l, ',') for l in lines[2:end] if !isempty(strip(l))], col
end

function main()
    rrows, rcol = read_csv(joinpath(DIR, "results.csv"))
    lrows, lcol = read_csv(joinpath(DIR, "landscape.csv"))

    key(f, c) = (f[c["family"]], f[c["n"]], f[c["instance"]], f[c["model"]])
    disp = Dict(key(f, lcol) => hypot(parse(Float64, f[lcol["dgamma"]]),
                                      parse(Float64, f[lcol["dbeta"]]))
                for f in lrows if f[lcol["tag"]] == "raw")

    groups = [
        ("analytical",        ["paper"],                          :black,      :diamond),
        ("Triangle, unfitted", ["tri_init"],                      :lightsalmon, :circle),
        ("Triangle, fitted",  ["tri_25", "tri_50", "tri_fit"],    :darkred,    :circle),
        ("Gaussian, unfitted", ["norm_init"],                      :lightskyblue, :utriangle),
        ("Gaussian, fitted",   ["norm_25", "norm_50", "norm_fit"], :navy,       :utriangle),
    ]
    bymodel = Dict(m => (Float64[], Float64[], Float64[])  # mse, disp, regret
                   for (_, ms, _, _) in groups for m in ms)
    for f in rrows
        mdl = f[rcol["model"]]
        haskey(bymodel, mdl) || continue
        mse, dsp, reg = bymodel[mdl]
        push!(mse, parse(Float64, f[rcol["mse_entry"]]))
        push!(dsp, disp[key(f, rcol)])
        push!(reg, parse(Float64, f[rcol["regret_raw"]]))
    end
    @assert all(length(v[1]) == 140 for v in values(bymodel)) "expected 140 rows per model"

    fig = Figure(size = (1150, 440), fontsize = 16)
    ax1 = Axis(fig[1, 1]; xscale = log10,
               xlabel = "fit objective: entrywise MSE vs empirical N",
               ylabel = "regret at p = 1",
               title = "fit quality does not predict regret")
    ax2 = Axis(fig[1, 2];
               xlabel = "argmax displacement ‖(Δγ, Δβ)‖ from empirical proxy",
               ylabel = "regret at p = 1",
               title = "argmax displacement does")
    for (label, models, color, marker) in groups
        mse = vcat((bymodel[m][1] for m in models)...)
        dsp = vcat((bymodel[m][2] for m in models)...)
        reg = vcat((bymodel[m][3] for m in models)...)
        scatter!(ax1, mse, reg; color, marker, markersize = 7, alpha = 0.55, label)
        scatter!(ax2, dsp, reg; color, marker, markersize = 7, alpha = 0.55, label)
    end
    axislegend(ax1; position = :lt, labelsize = 12)

    out = joinpath(DIR, "fitted_shape_paradox.png")
    save(out, fig)
    println("wrote $out")
    outpdf = joinpath(DIR, "fitted_shape_paradox.pdf")
    save(outpdf, fig)
    println("wrote $outpdf")
end

main()
