#=
E2.4 analysis: does entrywise-MSE fit quality predict parameter-setting regret,
and does the dynamics-weighted model error ‖Q₁_model − Q₁_emp‖ do better?

Reads results.csv (written by run.jl), prints a summary (also saved to
analysis_summary.txt) and writes figures:
  fig_regret_vs_mse.png      — regret vs the fit objective, raw and calibrated
  fig_regret_vs_eps.png      — regret vs dynamics-weighted error, raw and calibrated
  fig_fit_trajectories.png   — mean MSE and mean regret along the fit trajectory

Run from the repo root:
  julia --project research/experiments/012_fitted-shape-paradox/analyze.jl
=#

using Statistics: mean, median, std
using Printf: @sprintf
using CairoMakie

const DIR = @__DIR__
const RESULTS = joinpath(DIR, "results.csv")

# ---------------------------------------------------------------- CSV parsing

function read_results(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    rows = NamedTuple[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        f = split(line, ',')
        @assert length(f) == length(header) "ragged CSV row: $line"
        d = Dict(zip(header, f))
        push!(rows, (
            family = String(d["family"]),
            n = parse(Int, d["n"]),
            instance = parse(Int, d["instance"]),
            m = parse(Int, d["m"]),
            model = String(d["model"]),
            mse_entry = parse(Float64, d["mse_entry"]),
            eps_raw_mean = parse(Float64, d["eps_raw_mean"]),
            eps_raw_argmax = parse(Float64, d["eps_raw_argmax"]),
            eps_cal_mean = parse(Float64, d["eps_cal_mean"]),
            eps_cal_argmax = parse(Float64, d["eps_cal_argmax"]),
            ceil_p1 = parse(Float64, d["ceil_p1"]),
            ar_raw = parse(Float64, d["ar_raw"]),
            regret_raw = parse(Float64, d["regret_raw"]),
            ar_cal = parse(Float64, d["ar_cal"]),
            regret_cal = parse(Float64, d["regret_cal"]),
        ))
    end
    return rows
end

# ------------------------------------------------------------------- Spearman

"Average ranks with tie handling."
function tiedrank(v::AbstractVector{<:Real})
    p = sortperm(v)
    r = similar(v, Float64)
    i = 1
    while i <= length(v)
        j = i
        while j < length(v) && v[p[j+1]] == v[p[i]]
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

function spearman(x, y)
    @assert length(x) == length(y) && length(x) >= 3
    rx, ry = tiedrank(x), tiedrank(y)
    sx, sy = std(rx), std(ry)
    (sx == 0 || sy == 0) && return NaN  # constant ranks: correlation undefined
    return mean((rx .- mean(rx)) .* (ry .- mean(ry))) / (sx * sy)
end

# ------------------------------------------------------------------- analysis

const MODELS = ["paper",
                "tri_init", "tri_25", "tri_50", "tri_fit",
                "norm_init", "norm_25", "norm_50", "norm_fit"]

function main()
    rows = read_results(RESULTS)
    out = IOBuffer()
    prnt(s="") = (println(out, s); println(s))

    instances = unique((r.family, r.n, r.instance) for r in rows)
    byinst = Dict(k => Dict(r.model => r
                            for r in rows
                            if (r.family, r.n, r.instance) == k)
                  for k in instances)
    @assert all(length(v) == 10 for v in values(byinst)) "expected 10 models per instance"

    prnt("E2.4 analysis — $(length(instances)) instances, $(length(MODELS)) non-emp models each")
    prnt()

    # --- per-model summary table -------------------------------------------
    prnt("Per-model means over all instances:")
    prnt(rpad("model", 10) * join(rpad.(["mse_entry", "eps_raw_am", "eps_cal_am",
                                         "regret_raw", "regret_cal"], 12)))
    for mdl in ["emp"; MODELS]
        sel = [byinst[k][mdl] for k in instances]
        prnt(rpad(mdl, 10) * join(rpad.([
            @sprintf("%.3e", mean(r.mse_entry for r in sel)),
            @sprintf("%.4f", mean(r.eps_raw_argmax for r in sel)),
            @sprintf("%.4f", mean(r.eps_cal_argmax for r in sel)),
            @sprintf("%.4f", mean(r.regret_raw for r in sel)),
            @sprintf("%.4f", mean(r.regret_cal for r in sel))], 12)))
    end
    prnt()

    # --- per-instance Spearman across the model spectrum --------------------
    pairs = [
        ("regret_raw ~ mse_entry",      r -> r.mse_entry,      r -> r.regret_raw),
        ("regret_raw ~ eps_raw_mean",   r -> r.eps_raw_mean,   r -> r.regret_raw),
        ("regret_raw ~ eps_raw_argmax", r -> r.eps_raw_argmax, r -> r.regret_raw),
        ("regret_cal ~ mse_entry",      r -> r.mse_entry,      r -> r.regret_cal),
        ("regret_cal ~ eps_cal_mean",   r -> r.eps_cal_mean,   r -> r.regret_cal),
        ("regret_cal ~ eps_cal_argmax", r -> r.eps_cal_argmax, r -> r.regret_cal),
    ]
    prnt("Per-instance Spearman ρ across the 9-model spectrum (mean ± sd, median; pooled):")
    for (name, fx, fy) in pairs
        ρs = Float64[]
        for k in instances
            sel = [byinst[k][mdl] for mdl in MODELS]
            ρ = spearman(fx.(sel), fy.(sel))
            isnan(ρ) || push!(ρs, ρ)
        end
        pooled = spearman([fx(byinst[k][mdl]) for k in instances for mdl in MODELS],
                          [fy(byinst[k][mdl]) for k in instances for mdl in MODELS])
        prnt(rpad(name, 32) *
             @sprintf("%+.3f ± %.3f   med %+.3f   pooled %+.3f   (N=%d)",
                      mean(ρs), std(ρs), median(ρs), pooled, length(ρs)))
    end
    prnt()

    # --- H2: along each fit trajectory, MSE falls — does regret? ------------
    for (shape, traj) in [("Triangle", ["tri_init", "tri_25", "tri_50", "tri_fit"]),
                          ("Normal",   ["norm_init", "norm_25", "norm_50", "norm_fit"])]
        dmse  = Float64[]; draw = Float64[]; dcal = Float64[]
        worse_raw = 0; worse_cal = 0; mono = 0
        for k in instances
            sel = [byinst[k][mdl] for mdl in traj]
            mses = [r.mse_entry for r in sel]
            @assert issorted(mses, rev=true) || all(diff(mses) .<= 1e-15) "fit checkpoints must not increase MSE: $k $mses"
            push!(dmse, mses[end] - mses[1])
            push!(draw, sel[end].regret_raw - sel[1].regret_raw)
            push!(dcal, sel[end].regret_cal - sel[1].regret_cal)
            sel[end].regret_raw > sel[1].regret_raw + 1e-12 && (worse_raw += 1)
            sel[end].regret_cal > sel[1].regret_cal + 1e-12 && (worse_cal += 1)
            issorted([r.regret_raw for r in sel], rev=true) && (mono += 1)
        end
        Ninst = length(instances)
        prnt("$shape fit trajectory (init → fit), $Ninst instances:")
        prnt(@sprintf("  mean ΔMSE      = %+.3e (fit − init; negative = fit improved the objective)",
                      mean(dmse)))
        prnt(@sprintf("  mean Δregret_raw = %+.4f;  fit worse than init (raw): %d/%d;  regret monotone ↓: %d/%d",
                      mean(draw), worse_raw, Ninst, mono, Ninst))
        prnt(@sprintf("  mean Δregret_cal = %+.4f;  fit worse than init (cal): %d/%d",
                      mean(dcal), worse_cal, Ninst))
    end
    prnt()

    # --- calibration effect --------------------------------------------------
    for mdl in MODELS
        sel = [byinst[k][mdl] for k in instances]
        better = count(r.regret_cal < r.regret_raw - 1e-12 for r in sel)
        prnt(@sprintf("calibration: %-10s mean regret %.4f → %.4f (cal better on %d/%d)",
                      mdl, mean(r.regret_raw for r in sel),
                      mean(r.regret_cal for r in sel), better, length(sel)))
    end

    write(joinpath(DIR, "analysis_summary.txt"), String(take!(out)))

    # ------------------------------------------------------------- figures
    colors = Dict("paper" => :black,
                  "tri_init" => :lightsalmon, "tri_25" => :salmon,
                  "tri_50" => :tomato, "tri_fit" => :darkred,
                  "norm_init" => :lightskyblue, "norm_25" => :deepskyblue,
                  "norm_50" => :dodgerblue, "norm_fit" => :navy)

    function scatterfig(fname, xfun_raw, xfun_cal, xlabel)
        fig = Figure(size = (1100, 460))
        for (col, xf, yf, title) in [(1, xfun_raw, r -> r.regret_raw, "raw N"),
                                     (2, xfun_cal, r -> r.regret_cal, "slice-calibrated N")]
            ax = Axis(fig[1, col]; xscale = log10, xlabel, ylabel = "regret (p=1)",
                      title)
            for mdl in MODELS
                sel = [byinst[k][mdl] for k in instances]
                scatter!(ax, max.(xf.(sel), 1e-12), yf.(sel);
                         color = colors[mdl], markersize = 6, label = mdl)
            end
            col == 2 && axislegend(ax; position = :rt, labelsize = 9, rowgap = 0)
        end
        save(joinpath(DIR, fname), fig)
    end

    scatterfig("fig_regret_vs_mse.png",
               r -> r.mse_entry, r -> r.mse_entry,
               "entrywise MSE of fit (slice-normalized)")
    scatterfig("fig_regret_vs_eps.png",
               r -> r.eps_raw_argmax, r -> r.eps_cal_argmax,
               "dynamics-weighted error ‖Q₁_model − Q₁_emp‖ at emp argmax")

    # fit-trajectory figure: mean mse and mean regret per checkpoint
    fig = Figure(size = (1100, 420))
    ax1 = Axis(fig[1, 1]; ylabel = "mean entrywise MSE", yscale = log10,
               xticks = (1:4, ["init", "25%", "50%", "fit"]),
               title = "fit objective along trajectory")
    ax2 = Axis(fig[1, 2]; ylabel = "mean regret (p=1)",
               xticks = (1:4, ["init", "25%", "50%", "fit"]),
               title = "parameter-setting regret along trajectory")
    for (shape, traj, color) in [("Triangle", ["tri_init", "tri_25", "tri_50", "tri_fit"], :darkred),
                                 ("Normal", ["norm_init", "norm_25", "norm_50", "norm_fit"], :navy)]
        ms = [mean(byinst[k][mdl].mse_entry for k in instances) for mdl in traj]
        rr = [mean(byinst[k][mdl].regret_raw for k in instances) for mdl in traj]
        rc = [mean(byinst[k][mdl].regret_cal for k in instances) for mdl in traj]
        scatterlines!(ax1, 1:4, ms; color, label = shape)
        scatterlines!(ax2, 1:4, rr; color, label = "$shape raw")
        scatterlines!(ax2, 1:4, rc; color, linestyle = :dash, label = "$shape calibrated")
    end
    axislegend(ax1; position = :rt)
    axislegend(ax2; position = :rt, labelsize = 10)
    save(joinpath(DIR, "fig_fit_trajectories.png"), fig)

    println("\nfigures + analysis_summary.txt written to $DIR")
end

main()
