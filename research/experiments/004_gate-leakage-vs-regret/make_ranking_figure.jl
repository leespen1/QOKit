#=
Headline figure + regret table for §5.3 ("compression error and parameter
quality"). Synthesizes the exact-compression regret and the leakage→regret
ranking across the seven graph families, with the honest scope built in.

Data sources (read-only):
  004_gate-leakage-vs-regret/results.csv   n=12,14; leakage/overlap/regret at
                                            the PROXY-CHOSEN angles, p=1 and p=3
                                            (ar_emp = exact-compression proxy)
  010_scaleup-ranking/results_task*.csv     n=16,18; ar_exact = exact compression

Outputs (next to this script; copy ranking.png into the paper repo's
Figures/generated/ and \input regret_table.tex):
  ranking.png        two panels:
    (a) leakage (at proxy-chosen angles) vs regret at p=1, n=12 & 14, 7 families,
        regret ±SE, Spearman ρ annotated — the ranking holds.
    (b) Spearman ρ(leakage,regret) across regimes — strong at p=1, gone at p=3:
        the ranking is a p=1 statement, stated honestly.
  regret_table.tex   exact-compression regret (mean ±SE), family × n, p=1 & p=3.

NOTE on the two meanings of "leakage": here leakage is measured at the proxy's
own chosen angles (exp 004), where the proxy shrinks γ on dense graphs — so
ER(0.5) leaks LEAST and 3-regular MOST. This is the OPPOSITE family ordering
from §5.2's fixed-angle density law. Both are correct; they are different
measurements. The panel titles say "at proxy-chosen angles" to keep them apart.

Run from the repo root:
  julia --project research/experiments/004_gate-leakage-vs-regret/make_ranking_figure.jl
=#

using Statistics: mean, std
using Printf: @sprintf
using CairoMakie

const RES  = abspath(joinpath(@__DIR__, ".."))
const DIR  = @__DIR__
const FAMS = ["ER(0.5)", "ER(0.25)", "BA(k=2)", "BA(k=4)",
              "WS(k=4;b=0.1)", "WS(k=4;b=0.5)", "3-regular"]
const FAMLABEL = Dict("ER(0.5)"=>"ER(0.5)", "ER(0.25)"=>"ER(0.25)",
                      "BA(k=2)"=>"BA(2)", "BA(k=4)"=>"BA(4)",
                      "WS(k=4;b=0.1)"=>"WS(0.1)", "WS(k=4;b=0.5)"=>"WS(0.5)",
                      "3-regular"=>"3-reg")

function readcsv(path)
    lines = readlines(path)
    header = String.(split(lines[1], ','))
    cols = Dict{String,Vector{String}}(h => String[] for h in header)
    for ln in lines[2:end]
        isempty(ln) && continue
        f = split(ln, ',')
        for (j, h) in enumerate(header); push!(cols[h], String(f[j])); end
    end
    return cols
end
fcol(cols, name) = parse.(Float64, cols[name])
icol(cols, name) = parse.(Int, cols[name])

function load010()
    cols = readcsv(joinpath(RES, "010_scaleup-ranking", "results_task1.csv"))
    cols = Dict(k => copy(v) for (k, v) in cols)
    for t in 2:14
        ct = readcsv(joinpath(RES, "010_scaleup-ranking", "results_task$t.csv"))
        for k in keys(cols); append!(cols[k], ct[k]); end
    end
    return cols
end

cellidx(cols, fam, n) =
    findall(i -> cols["family"][i] == fam && parse(Int, cols["n"][i]) == n,
            eachindex(cols["family"]))

function spearman(x, y)
    rx = Float64.(sortperm(sortperm(x)))
    ry = Float64.(sortperm(sortperm(y)))
    mx, my = mean(rx), mean(ry)
    sum((rx .- mx) .* (ry .- my)) /
        sqrt(sum((rx .- mx) .^ 2) * sum((ry .- my) .^ 2))
end

mse(v) = (mean(v), std(v) / sqrt(length(v)))

# (mean leakage, mean regret, regret SE) per family for a given dataset/regime.
function family_stats(cols, n, leakcol, ceilcol, arcol)
    leak = Float64[]; reg = Float64[]; se = Float64[]
    for fam in FAMS
        idx = cellidx(cols, fam, n)
        @assert !isempty(idx) "no rows for $fam n=$n"
        push!(leak, mean(fcol(cols, leakcol)[idx]))
        r = fcol(cols, ceilcol)[idx] .- fcol(cols, arcol)[idx]
        m, s = mse(r); push!(reg, m); push!(se, s)
    end
    return leak, reg, se
end

function main()
    c4  = readcsv(joinpath(DIR, "results.csv"))
    c10 = load010()

    # --- fail-fast sanity ---
    @assert sort(unique(c4["family"])) == sort(FAMS) "exp004 families differ"
    @assert sort(unique(icol(c4, "n"))) == [12, 14]
    @assert sort(unique(icol(c10, "n"))) == [16, 18]
    for fam in FAMS, n in (12, 14)
        @assert length(cellidx(c4, fam, n)) == 30 "exp004 $fam n=$n not 30 instances"
    end
    for fam in FAMS, n in (16, 18)
        @assert length(cellidx(c10, fam, n)) == 20 "exp010 $fam n=$n not 20 instances"
    end

    colors = Makie.wong_colors()[1:7]

    fig = Figure(size = (1150, 470), fontsize = 16)

    # ---- panel (a): p=1 leakage vs regret, n=12 & 14 ----
    axa = Axis(fig[1, 1];
        xlabel = "mean leakage Σλ at proxy-chosen angles",
        ylabel = "mean regret  (ceiling − proxy AR)",
        title = "(a) p = 1: leakage ranks regret")
    rho_legend = String[]
    for (n, marker) in [(12, :circle), (14, :rect)]
        leak, reg, se = family_stats(c4, n, "sum_leakage_p1", "ceil_p1", "ar_emp_p1")
        errorbars!(axa, leak, reg, se; color = :gray70, whiskerwidth = 6)
        scatter!(axa, leak, reg; color = colors, marker = marker, markersize = 14,
                 strokecolor = :black, strokewidth = 0.5)
        push!(rho_legend, "n=$n: ρ=$(round(spearman(leak, reg), digits = 2))")
    end
    text!(axa, 0.02, 0.97; text = join(rho_legend, "\n"), space = :relative,
          align = (:left, :top), fontsize = 14)
    # family color legend + marker legend
    famelems = [MarkerElement(color = colors[i], marker = :circle, markersize = 12)
                for i in 1:7]
    nelems = [MarkerElement(color = :black, marker = m, markersize = 12)
              for m in (:circle, :rect)]
    Legend(fig[1, 1], [famelems, nelems],
           [[FAMLABEL[f] for f in FAMS], ["n=12", "n=14"]],
           ["family", "size"];
           tellwidth = false, tellheight = false, halign = :right, valign = :bottom,
           labelsize = 11, titlesize = 12, patchsize = (12, 12), rowgap = 0,
           margin = (6, 6, 6, 6), framevisible = true, nbanks = 1)

    # ---- panel (b): the fade — ρ across regimes ----
    axb = Axis(fig[1, 2];
        ylabel = "Spearman ρ(leakage, regret)  over 7 families",
        title = "(b) the ranking is a p = 1 statement",
        xticks = (1:6, ["p1\nn=12", "p1\nn=14", "p3\nn=12", "p3\nn=14",
                        "p3\nn=16", "p3\nn=18"]))
    rhos = Float64[]
    # p=1 and p=3 at n=12,14 from exp004; p=3 at n=16,18 from exp010
    for n in (12, 14)
        leak, reg, _ = family_stats(c4, n, "sum_leakage_p1", "ceil_p1", "ar_emp_p1")
        push!(rhos, spearman(leak, reg))
    end
    for n in (12, 14)
        leak, reg, _ = family_stats(c4, n, "sum_leakage_p3", "ceil_p3", "ar_emp_p3")
        push!(rhos, spearman(leak, reg))
    end
    for n in (16, 18)
        leak, reg, _ = family_stats(c10, n, "sum_leakage_p3", "ceil_p3", "ar_exact_p3")
        push!(rhos, spearman(leak, reg))
    end
    barcolors = [:seagreen, :seagreen, :firebrick, :firebrick, :firebrick, :firebrick]
    barplot!(axb, 1:6, rhos; color = barcolors)
    hlines!(axb, [0.0]; color = :black, linewidth = 0.8)
    ylims!(axb, -0.2, 1.05)
    for (i, r) in enumerate(rhos)
        text!(axb, i, r + (r ≥ 0 ? 0.03 : -0.03); text = string(round(r, digits = 2)),
              align = (:center, r ≥ 0 ? :bottom : :top), fontsize = 12)
    end

    out = joinpath(DIR, "ranking.png")
    save(out, fig)
    println("wrote $out")
    outpdf = joinpath(DIR, "ranking.pdf")
    save(outpdf, fig)
    println("wrote $outpdf")

    # ---- regret table (LaTeX fragment) ----
    fmt(x) = @sprintf("%.3f", x)
    function regret_row(fam)
        v1 = String[]; v3 = String[]
        for n in (12, 14)
            idx = cellidx(c4, fam, n)
            push!(v1, fmt(mean(fcol(c4, "ceil_p1")[idx] .- fcol(c4, "ar_emp_p1")[idx])))
            push!(v3, fmt(mean(fcol(c4, "ceil_p3")[idx] .- fcol(c4, "ar_emp_p3")[idx])))
        end
        for n in (16, 18)
            idx = cellidx(c10, fam, n)
            push!(v1, fmt(mean(fcol(c10, "ceil_p1")[idx] .- fcol(c10, "ar_exact_p1")[idx])))
            push!(v3, fmt(mean(fcol(c10, "ceil_p3")[idx] .- fcol(c10, "ar_exact_p3")[idx])))
        end
        return v1, v3
    end
    open(joinpath(DIR, "regret_table.tex"), "w") do io
        println(io, "% Auto-generated by make_ranking_figure.jl (exp 004/010). Do not hand-edit.")
        println(io, "% Exact-compression regret (mean over instances; per-cell SE ~0.001).")
        println(io, "\\begin{tabular}{lcccc|cccc}")
        println(io, "  \\toprule")
        println(io, "  & \\multicolumn{4}{c|}{regret, \$p=1\$} & \\multicolumn{4}{c}{regret, \$p=3\$} \\\\")
        println(io, "  family & \$n{=}12\$ & \$14\$ & \$16\$ & \$18\$ & \$n{=}12\$ & \$14\$ & \$16\$ & \$18\$ \\\\")
        println(io, "  \\midrule")
        for fam in FAMS
            v1, v3 = regret_row(fam)
            println(io, "  " * FAMLABEL[fam] * " & " * join(v1, " & ") * " & " *
                        join(v3, " & ") * " \\\\")
        end
        println(io, "  \\bottomrule")
        println(io, "\\end{tabular}")
    end
    println("wrote ", joinpath(DIR, "regret_table.tex"))

    # ---- console summary + fail-fast on the headline number ----
    println("\nSpearman ρ by regime: ", round.(rhos, digits = 3))
    @assert isapprox(rhos[1], 0.96; atol = 0.03) "n=12 p=1 ρ drifted from 0.96: $(rhos[1])"
    @assert rhos[3] < 0.3 && rhos[4] < 0.3 "p=3 n=12/14 ρ should be ~0 (the fade)"
    println("OK: headline ρ reproduced; p=3 fade confirmed.")
end

main()
