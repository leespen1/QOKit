#=
Analysis for E014. Reads results.csv and asks, at p=1 and p=3 separately:
what governs exact-compression regret — fidelity (1−overlap), argmax
displacement, or true-landscape robustness (flat-peak fraction)?

Framing question (STATUS): at p=3 proxy-chosen schedules equalize leakage/fidelity
across families, so if regret there tracks landscape robustness rather than
fidelity, depth regret is an argmax/landscape phenomenon, not a fidelity one.

Run: julia --project research/experiments/014_argmax-robustness/analyze.jl
=#
using Statistics: mean, std, cor

const DIR = @__DIR__
function readcsv(path)
    lines = readlines(path); h = String.(split(lines[1], ','))
    cols = Dict(c => String[] for c in h)
    for ln in lines[2:end]; isempty(ln) && continue; f = split(ln, ',')
        for (j, c) in enumerate(h); push!(cols[c], String(f[j])); end; end
    cols
end
fc(c, k) = parse.(Float64, c[k])
function spearman(x, y)
    rx = Float64.(sortperm(sortperm(x))); ry = Float64.(sortperm(sortperm(y))); cor(rx, ry)
end
cv(v) = std(v) / mean(v)   # coefficient of variation

c = readcsv(joinpath(DIR, "results.csv"))
fam = c["family"]; P = parse.(Int, c["p"])
reg = fc(c, "regret"); ov = fc(c, "overlap"); disp = fc(c, "argmax_disp"); rob = fc(c, "robust_frac")
FAMS = unique(fam)

for p in (1, 3)
    s = P .== p
    fdef = 1 .- ov[s]                       # fidelity deficit
    println("\n================  p = $p  ================")
    println("POOLED Spearman(regret, ·) over $(count(s)) instances:")
    println("   fidelity deficit (1−overlap) : ", round(spearman(reg[s], fdef), digits=3))
    println("   argmax displacement          : ", round(spearman(reg[s], disp[s]), digits=3))
    println("   landscape robustness (flat%) : ", round(spearman(reg[s], rob[s]), digits=3), "  (expect negative)")
    # family means (7 points)
    fmean(v) = [mean(v[s .& (fam .== f)]) for f in FAMS]
    mreg = fmean(reg); mfdef = [mean((1 .- ov[s .& (fam .== f)])) for f in FAMS]
    mdisp = fmean(disp); mrob = fmean(rob); mov = fmean(ov)
    println("FAMILY-MEAN Spearman over 7 families:")
    println("   fidelity deficit : ", round(spearman(mreg, mfdef), digits=3))
    println("   argmax disp      : ", round(spearman(mreg, mdisp), digits=3))
    println("   robustness       : ", round(spearman(mreg, mrob), digits=3))
    println("Cross-family spread (CV) of each predictor — is fidelity 'equalized' here?")
    println("   CV(overlap)=", round(cv(mov), digits=3),
            "  CV(robust)=", round(cv(mrob), digits=3),
            "  CV(regret)=", round(cv(mreg), digits=3))
    println("per-family means:")
    for (i, f) in enumerate(FAMS)
        println("   ", rpad(f, 14), " regret=", round(mreg[i], digits=4),
                "  overlap=", round(mov[i], digits=4),
                "  robust=", round(mrob[i], digits=4),
                "  disp=", round(mdisp[i], digits=4))
    end
end
