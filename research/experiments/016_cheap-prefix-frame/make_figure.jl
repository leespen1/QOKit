#=
Figure for E016. Per ramp regime:
  (a) captured energy of the full P=20 trajectory by three frames — cost-class
      (dim m+1, 0 layers), optimal top-4 PCA (oracle), cheap top-4 PCA of the
      first 5 layers. The oracle shows the headroom; the cheap prefix falls below
      even the zero-cost cost-class frame at working ramps.
  (b) largest principal angle between the cheap-prefix and full top-4 subspaces —
      near-orthogonal (71–88°) and growing with ramp: the subspace rotates over
      depth, so an early prefix cannot discover the late frame.

Writes prefix_frame.png next to this script.
Run: julia --project research/experiments/016_cheap-prefix-frame/make_figure.jl
=#
using Statistics: mean
using CairoMakie

const DIR = @__DIR__
function readcsv(path)
    lines = readlines(path); h = String.(split(lines[1], ','))
    cols = Dict(c => String[] for c in h)
    for ln in lines[2:end]; isempty(ln) && continue; f = split(ln, ',')
        for (j, c) in enumerate(h); push!(cols[c], String(f[j])); end; end
    cols
end
fc(c, k) = parse.(Float64, c[k])

function main()
    c = readcsv(joinpath(DIR, "results.csv"))
    ramp = c["ramp"]; Ecc = fc(c,"E_cc"); Eful = fc(c,"E_full_d"); Epre = fc(c,"E_prefix_d")
    ang = fc(c,"max_angle_deg")
    ramps = ["small", "moderate", "large"]
    mean_by(v) = [mean(v[ramp .== r]) for r in ramps]

    fig = Figure(size = (1080, 440), fontsize = 15)
    axa = Axis(fig[1, 1]; ylabel = "captured energy of full p=20 trajectory",
        title = "(a) the headroom exists, but the cheap prefix misses it",
        xticks = (1:3, ramps))
    ccm, fum, prm = mean_by(Ecc), mean_by(Eful), mean_by(Epre)
    xs = repeat(1:3, 3)
    grp = vcat(fill(1, 3), fill(2, 3), fill(3, 3))
    ys = vcat(fum, ccm, prm)        # oracle, cost-class, prefix
    cols3 = [:seagreen, :slategray, :firebrick]
    barplot!(axa, xs, ys; dodge = grp, color = cols3[grp])
    ylims!(axa, 0.6, 1.02)
    elems = [PolyElement(color = cols3[i]) for i in 1:3]
    Legend(fig[1, 1], elems, ["optimal PCA d=4 (oracle, needs all p)",
        "cost-class (dim m+1, 0 layers)", "cheap prefix d=4 (first 5 layers)"];
        tellwidth = false, tellheight = false, halign = :left, valign = :bottom,
        labelsize = 11, framevisible = true, margin = (8,8,8,8))

    axb = Axis(fig[1, 2]; ylabel = "largest principal angle  ∠(prefix, full)  [deg]",
        title = "(b) the 4-dim subspace rotates over depth",
        xticks = (1:3, ramps))
    barplot!(axb, 1:3, mean_by(ang); color = :purple)
    hlines!(axb, [90.0]; color = :black, linestyle = :dash)
    text!(axb, 0.5, 90.5; text = "90° = orthogonal", align = (:left, :bottom), fontsize = 12)
    ylims!(axb, 0, 100)

    out = joinpath(DIR, "prefix_frame.png"); save(out, fig); println("wrote $out")
end
main()
