# E017: statistics hardening for the paper's headline correlations.
# Analysis-only; consumes committed CSVs from E004/E010/E011/E012/E014.
# Outputs results.csv (one row per statistic) and prints a summary.
#
# Reproduce: julia --project research/experiments/017_statistics-hardening/run.jl
using Random

const HERE = @__DIR__
const EXP = dirname(HERE)

# ---------- tiny CSV reader (schemas here have no quoted commas) ----------
function readcsv(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    cols = Dict(String(h) => Vector{String}() for h in header)
    for ln in lines[2:end]
        isempty(strip(ln)) && continue
        parts = split(ln, ',')
        length(parts) == length(header) ||
            error("bad row in $path: $(length(parts)) fields vs $(length(header))")
        for (h, v) in zip(header, parts)
            push!(cols[String(h)], String(v))
        end
    end
    return cols
end
fl(col) = parse.(Float64, col)

# ---------- tie-averaged Spearman ----------
function ranks(x::Vector{Float64})
    p = sortperm(x)
    r = similar(x)
    i = 1
    while i <= length(x)
        j = i
        while j < length(x) && x[p[j+1]] == x[p[i]]
            j += 1
        end
        avg = (i + j) / 2
        for k in i:j
            r[p[k]] = avg
        end
        i = j + 1
    end
    return r
end
function pearson(x, y)
    mx, my = sum(x)/length(x), sum(y)/length(y)
    sx = sqrt(sum(abs2, x .- mx)); sy = sqrt(sum(abs2, y .- my))
    (sx == 0 || sy == 0) && return NaN
    return sum((x .- mx) .* (y .- my)) / (sx * sy)
end
spearman(x, y) = pearson(ranks(collect(Float64, x)), ranks(collect(Float64, y)))

# ---------- bootstrap helpers ----------
"Percentile CI from a statistic over resamples."
function ci(vals; lo=0.025, hi=0.975)
    v = sort(filter(!isnan, vals))
    isempty(v) && return (NaN, NaN)
    return (v[max(1, ceil(Int, lo*length(v)))], v[min(length(v), ceil(Int, hi*length(v)))])
end

"Instance-level bootstrap CI of Spearman(x,y)."
function boot_rho(rng, x, y; B=5000)
    n = length(x)
    vals = Vector{Float64}(undef, B)
    for b in 1:B
        idx = rand(rng, 1:n, n)
        vals[b] = spearman(x[idx], y[idx])
    end
    return ci(vals)
end

"Cluster bootstrap: resample whole cells (vectors of indices) with replacement."
function boot_rho_cluster(rng, x, y, cells; B=5000)
    vals = Vector{Float64}(undef, B)
    for b in 1:B
        idx = Int[]
        for _ in 1:length(cells)
            append!(idx, cells[rand(rng, 1:length(cells))])
        end
        vals[b] = spearman(x[idx], y[idx])
    end
    return ci(vals)
end

"Mean within-cell Spearman +- SE, and family-demeaned pooled Spearman."
function within_cell(x, y, cellid)
    ids = unique(cellid)
    rhos = Float64[]
    xd = copy(x); yd = copy(y)
    for c in ids
        sel = findall(==(c), cellid)
        length(sel) >= 4 || continue
        r = spearman(x[sel], y[sel])
        isnan(r) || push!(rhos, r)  # NaN = a zero-variance cell; skip it
        xd[sel] .-= sum(x[sel])/length(sel)
        yd[sel] .-= sum(y[sel])/length(sel)
    end
    m = sum(rhos)/length(rhos)
    se = sqrt(sum(abs2, rhos .- m) / (length(rhos)-1) / length(rhos))
    return m, se, spearman(xd, yd), length(rhos)
end

results = String["statistic,value,ci_lo,ci_hi,detail"]
row(name, val, lo, hi, detail) =
    push!(results, "$name,$(round(val; digits=4)),$(round(lo; digits=4)),$(round(hi; digits=4)),$detail")

rng = MersenneTwister(20260722)

# ================= E014: regret predictors at p=1 and p=3 =================
e14 = readcsv(joinpath(EXP, "014_argmax-robustness", "results.csv"))
for p in ("1", "3")
    sel = findall(==(p), e14["p"])
    regret = fl(e14["regret"][sel])
    cellid = [e14["family"][i] * "|" * e14["n"][i] for i in sel]
    cells_map = Dict{String,Vector{Int}}()
    for (k, c) in enumerate(cellid)
        push!(get!(cells_map, c, Int[]), k)
    end
    cells = collect(values(cells_map))
    for (label, vals) in (
        ("argmax_disp", fl(e14["argmax_disp"][sel])),
        ("fidelity_deficit", 1 .- fl(e14["overlap"][sel])),
        ("robust_frac", fl(e14["robust_frac"][sel])),
    )
        rho = spearman(regret, vals)
        lo, hi = boot_rho(rng, regret, vals)
        clo, chi = boot_rho_cluster(rng, regret, vals, cells)
        wm, wse, demeaned, ncells = within_cell(regret, vals, cellid)
        row("e014_p$(p)_rho_$(label)", rho, lo, hi, "pooled n=$(length(sel)) instance-bootstrap")
        row("e014_p$(p)_rho_$(label)_cluster", rho, clo, chi, "cluster bootstrap over $(length(cells)) family-n cells")
        row("e014_p$(p)_rho_$(label)_withincell", wm, wm - 1.96*wse, wm + 1.96*wse, "mean of $(ncells) within-cell rhos (10 inst each)")
        row("e014_p$(p)_rho_$(label)_demeaned", demeaned, NaN, NaN, "cell-demeaned pooled Spearman")
    end
end

# ================= E012: mismatch norms vs regret =================
e12 = readcsv(joinpath(EXP, "012_fitted-shape-paradox", "results.csv"))
sel = findall(!=("emp"), e12["model"])
regret = fl(e12["regret_raw"][sel])
cellid = [e12["family"][i] * "|" * e12["n"][i] for i in sel]
for (label, vals) in (
    ("mse_entry", fl(e12["mse_entry"][sel])),
    ("eps_raw_argmax", fl(e12["eps_raw_argmax"][sel])),
)
    rho = spearman(regret, vals)
    lo, hi = boot_rho(rng, regret, vals)
    wm, wse, demeaned, ncells = within_cell(regret, vals, cellid)
    row("e012_rho_$(label)", rho, lo, hi, "pooled over $(length(sel)) (instance,model) rows")
    row("e012_rho_$(label)_withincell", wm, wm - 1.96*wse, wm + 1.96*wse, "mean of $(ncells) within-cell rhos")
    row("e012_rho_$(label)_demeaned", demeaned, NaN, NaN, "cell-demeaned pooled Spearman")
end

# ================= E004: 7-point family ranking with bootstrap CI ==========
e04 = readcsv(joinpath(EXP, "004_gate-leakage-vs-regret", "results.csv"))
for nval in ("12", "14")
    local sel = findall(==(nval), e04["n"])
    fams = unique(e04["family"][sel])
    length(fams) == 7 || error("expected 7 families at n=$nval, got $(length(fams))")
    leak = fl(e04["sum_leakage_p1"][sel])
    reg = fl(e04["ceil_p1"][sel]) .- fl(e04["ar_emp_p1"][sel])
    fam = e04["family"][sel]
    fmean(v, f) = [sum(v[fam .== g]) / count(==(g), fam) for g in f]
    rho = spearman(fmean(leak, fams), fmean(reg, fams))
    B = 5000
    vals = Vector{Float64}(undef, B)
    for b in 1:B
        lb = Float64[]; rb = Float64[]
        for g in fams
            gi = findall(==(g), fam)
            idx = rand(rng, gi, length(gi))
            push!(lb, sum(leak[idx])/length(idx))
            push!(rb, sum(reg[idx])/length(idx))
        end
        vals[b] = spearman(lb, rb)
    end
    lo, hi = ci(vals)
    row("e004_p1_ranking_rho_n$(nval)", rho, lo, hi, "7 family means, instance bootstrap within family")
end

# ================= E011: depth-accumulation ratio spread ===================
# Paper method: per-instance Sum-lambda of the FULL p=30 run (E011) divided by
# the FULL p=20 run (E003) at identical family/instance/ramp endpoints, n=16.
function final_cum(dir, files, nval, pmax)
    out = Dict{String,Float64}()
    for f in files
        d = readcsv(joinpath(EXP, dir, f))
        for i in eachindex(d["family"])
            d["n"][i] == nval || continue
            d["ramp"][i] == "small" || continue
            parse(Int, d["layer"][i]) == pmax || continue
            out[d["family"][i] * "|" * d["instance"][i]] = parse(Float64, d["cum_leakage"][i])
        end
    end
    return out
end
cum30 = final_cum("011_depth-scaling",
    filter(x -> startswith(x, "results_task"), readdir(joinpath(EXP, "011_depth-scaling"))), "16", 30)
cum20 = final_cum("003_leakage-vs-overlap", ["results.csv"], "16", 20)
common = intersect(keys(cum30), keys(cum20))
isempty(common) && error("no matched (family, instance) runs between E011 and E003 at n=16")
ratios = [cum30[k] / cum20[k] for k in common]
m = sum(ratios)/length(ratios)
sd = sqrt(sum(abs2, ratios .- m)/(length(ratios)-1))
row("e011_ratio_p30_p20_small", m, m - sd, m + sd,
    "per-instance ratio over $(length(ratios)) matched n=16 small-ramp runs (CI cols = +-1 SD)")

# ================= E010: sampled-vs-exact regret sign test =================
rows10 = Dict{String,Vector{Float64}}()  # cell -> [sum_exact, sum_samp, count]
for f in filter(x -> startswith(x, "results_task"), readdir(joinpath(EXP, "010_scaleup-ranking")))
    d = readcsv(joinpath(EXP, "010_scaleup-ranking", f))
    for i in eachindex(d["family"])
        for p in ("p1", "p3")
            cell = join((d["family"][i], d["n"][i], p), "|")
            v = get!(rows10, cell, [0.0, 0.0, 0.0])
            ceilv = parse(Float64, d["ceil_$p"][i])
            v[1] += ceilv - parse(Float64, d["ar_exact_$p"][i])
            v[2] += ceilv - parse(Float64, d["ar_samp_$p"][i])
            v[3] += 1
        end
    end
end
lower = count(v -> v[2] < v[1], values(rows10))
total = length(rows10)
# one-sided binomial P(X >= lower | total, 1/2)
logbinom(n, k) = sum(log.(k+1:n)) - sum(log.(1:n-k))
pval = sum(exp(logbinom(total, k) + total*log(0.5)) for k in lower:total)
push!(results, "e010_sampled_lower_cells,$lower,$total,$(pval),one-sided binomial p in ci_hi col; cells with lower sampled regret out of total")

open(joinpath(HERE, "results.csv"), "w") do io
    for r in results
        println(io, r)
    end
end
println("Wrote $(length(results)-1) statistics to results.csv")
for r in results
    println(r)
end
