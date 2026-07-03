#=
Figure for E014. Two scatter panels over all instances, p=1 (blue) and p=3 (red):
  (a) regret vs fidelity deficit (1−overlap): a slope at p=1 that flattens to a
      cloud at p=3 (ρ: 0.39 → −0.02) — fidelity STOPS predicting regret at depth.
  (b) regret vs argmax displacement: a positive trend at BOTH depths
      (ρ: 0.76 → 0.63) — parameter-space distance keeps predicting regret.
Message: parameter-setting regret is an argmax-transfer quantity, decoupled from
state fidelity at depth.

Writes argmax_vs_fidelity.png next to this script.
Run: julia --project research/experiments/014_argmax-robustness/make_figure.jl
=#
using CairoMakie

const DIR = @__DIR__
function readcsv(path)
    lines = readlines(path); h = String.(split(lines[1], ','))
    cols = Dict(c => String[] for c in h)
    for ln in lines[2:end]; isempty(ln) && continue; f = split(ln, ',')
        for (j, cc) in enumerate(h); push!(cols[cc], String(f[j])); end; end
    cols
end
fcol(c, k) = parse.(Float64, c[k])
function spearman(x, y)
    rx = Float64.(sortperm(sortperm(x))); ry = Float64.(sortperm(sortperm(y)))
    mx, my = sum(rx)/length(rx), sum(ry)/length(ry)
    sum((rx.-mx).*(ry.-my)) / sqrt(sum((rx.-mx).^2)*sum((ry.-my).^2))
end

function main()
    c = readcsv(joinpath(DIR, "results.csv"))
    p = parse.(Int, c["p"]); reg = fcol(c, "regret")
    fdef = 1 .- fcol(c, "overlap"); disp = fcol(c, "argmax_disp")
    s1 = p .== 1; s3 = p .== 3

    fig = Figure(size = (1100, 470), fontsize = 15)
    axa = Axis(fig[1, 1]; xlabel = "fidelity deficit  1 − |⟨ψ|φ⟩|²",
        ylabel = "regret  (ceiling − proxy AR)",
        title = "(a) fidelity stops predicting regret at depth")
    scatter!(axa, fdef[s1], reg[s1]; color = (:navy, 0.55), markersize = 8, label = "p=1")
    scatter!(axa, fdef[s3], reg[s3]; color = (:firebrick, 0.55), markersize = 8, label = "p=3")
    text!(axa, 0.03, 0.97; space = :relative, align = (:left, :top), fontsize = 14,
        text = "Spearman ρ(regret, deficit)\n  p=1: $(round(spearman(reg[s1],fdef[s1]),digits=2))\n  p=3: $(round(spearman(reg[s3],fdef[s3]),digits=2))  ← decoupled")
    axislegend(axa; position = :rb)

    axb = Axis(fig[1, 2]; xlabel = "argmax displacement  (normalized angle distance)",
        ylabel = "regret  (ceiling − proxy AR)",
        title = "(b) argmax displacement predicts regret at both depths")
    scatter!(axb, disp[s1], reg[s1]; color = (:navy, 0.55), markersize = 8, label = "p=1")
    scatter!(axb, disp[s3], reg[s3]; color = (:firebrick, 0.55), markersize = 8, label = "p=3")
    text!(axb, 0.03, 0.97; space = :relative, align = (:left, :top), fontsize = 14,
        text = "Spearman ρ(regret, disp)\n  p=1: $(round(spearman(reg[s1],disp[s1]),digits=2))\n  p=3: $(round(spearman(reg[s3],disp[s3]),digits=2))")
    axislegend(axb; position = :rb)

    out = joinpath(DIR, "argmax_vs_fidelity.png")
    save(out, fig); println("wrote $out")
end
main()
