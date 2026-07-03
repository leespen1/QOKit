#=
E014 analysis — reads results_task*.csv and prints the summary tables that
back the README's Answer. Pure stdlib (plain CSV split; no quoted fields).

Run: julia research/experiments/014_fitted-shape-paradox/analyze.jl
=#

using Statistics

struct Row
    family::String
    n::Int
    instance::Int
    m::Int
    ceil::Float64
    variant::String
    mse::Float64
    ew_avg::Float64
    ew_star_raw::Float64
    ew_star_nrm::Float64
    ar_raw::Float64
    regret_raw::Float64
    ar_nrm::Float64
    regret_nrm::Float64
    pred::Float64
    normw::Float64
end

function load_rows(dir)
    rows = Row[]
    files = sort(filter(f -> occursin(r"^results_task\d+\.csv$", f), readdir(dir)))
    isempty(files) && error("no results_task*.csv in $dir")
    for f in files
        for (i, line) in enumerate(eachline(joinpath(dir, f)))
            i == 1 && continue
            p = split(line, ",")
            push!(rows, Row(p[1], parse(Int, p[2]), parse(Int, p[3]),
                            parse(Int, p[5]), parse(Float64, p[7]), p[8],
                            parse(Float64, p[9]), parse(Float64, p[10]),
                            parse(Float64, p[11]), parse(Float64, p[12]),
                            parse(Float64, p[13]), parse(Float64, p[14]),
                            parse(Float64, p[17]), parse(Float64, p[18]),
                            parse(Float64, p[21]), parse(Float64, p[22])))
        end
    end
    println("loaded $(length(rows)) rows from $(length(files)) task files\n")
    rows
end

function midranks(x)
    p = sortperm(x)
    r = zeros(length(x))
    i = 1
    while i <= length(x)
        j = i
        while j < length(x) && x[p[j + 1]] == x[p[i]]
            j += 1
        end
        for k in i:j
            r[p[k]] = (i + j) / 2
        end
        i = j + 1
    end
    r
end

spearman(x, y) = cor(midranks(x), midranks(y))

fmt(x; d=4) = isnan(x) ? "  nan " : lpad(string(round(x; digits=d)), 7)

function main()
    rows = load_rows(@__DIR__)
    cells = sort(unique((r.family, r.n) for r in rows))
    perts = filter(r -> startswith(r.variant, "pert:"), rows)

    # sanity: Thm-1 contractivity for exact N
    exact = filter(r -> r.variant == "exactN", rows)
    println("== Sanity ==")
    println("exactN max normw_raw = $(maximum(r.normw for r in exact)) (Thm 1 ⇒ ≤ 1)")
    println("exactN mean regret: raw $(fmt(mean(r.regret_raw for r in exact))), " *
            "nrm $(fmt(mean(r.regret_nrm for r in exact)))\n")

    # ---- Part A: matched-MSE perturbations --------------------------------
    println("== Part A: regret at matched entrywise MSE (pooled over cells) ==")
    println(rpad("shape", 9), "ε      mean mse_rel  mean ew_avg  mean regret_nrm  mean regret_raw")
    for shape in ["lowd", "highd", "nulled", "random"], ε in [0.01, 0.05, 0.2, 0.5]
        sel = filter(r -> r.variant == "pert:$shape:$ε", perts)
        isempty(sel) && continue
        println(rpad(shape, 9), rpad(ε, 7),
                fmt(mean(r.mse for r in sel)), "      ", fmt(mean(r.ew_avg for r in sel)),
                "     ", fmt(mean(r.regret_nrm for r in sel)),
                "         ", fmt(mean(r.regret_raw for r in sel)))
    end

    println("\nSpearman(regret_nrm, ·) over pert rows (+exactN anchor), per cell:")
    println(rpad("cell", 22), "ρ(ew_avg)  ρ(ew_star_nrm)  ρ(mse_rel)   [n rows]")
    pooled = Row[]
    for (fam, n) in cells
        sel = filter(r -> (startswith(r.variant, "pert:") || r.variant == "exactN") &&
                          r.family == fam && r.n == n, rows)
        append!(pooled, sel)
        reg = [r.regret_nrm for r in sel]
        println(rpad("$fam n=$n", 22),
                fmt(spearman(reg, [r.ew_avg for r in sel]); d=3), "    ",
                fmt(spearman(reg, [r.ew_star_nrm for r in sel]); d=3), "        ",
                fmt(spearman(reg, [r.mse for r in sel]); d=3), "   [", length(sel), "]")
    end
    reg = [r.regret_nrm for r in pooled]
    println(rpad("POOLED", 22),
            fmt(spearman(reg, [r.ew_avg for r in pooled]); d=3), "    ",
            fmt(spearman(reg, [r.ew_star_nrm for r in pooled]); d=3), "        ",
            fmt(spearman(reg, [r.mse for r in pooled]); d=3))

    # ---- Part B: fitted shapes vs PaperProxy ------------------------------
    println("\n== Part B: models (pooled means) ==")
    models = ["fit:triangle:mse", "fit:normal:mse", "fit:triangle:wgt",
              "fit:normal:wgt", "paper:empP", "paper:binP"]
    println(rpad("model", 18), "mse_rel  ew_avg   regret_nrm  regret_raw  normw(med)")
    for v in models
        sel = filter(r -> r.variant == v, rows)
        isempty(sel) && continue
        println(rpad(v, 18), fmt(mean(r.mse for r in sel)), "  ", fmt(mean(r.ew_avg for r in sel)),
                "  ", fmt(mean(r.regret_nrm for r in sel)), "     ", fmt(mean(r.regret_raw for r in sel)),
                "     ", fmt(median(r.normw for r in sel); d=2))
    end

    println("\nPairwise ranking concordance with regret_nrm across instances")
    println("(fraction of instance-pairs where the metric orders the two models")
    println(" the same way regret does; ties in regret skipped):")
    bykey = Dict{Tuple{String, Int, Int, String}, Row}()
    for r in rows
        bykey[(r.family, r.n, r.instance, r.variant)] = r
    end
    instances = sort(unique((r.family, r.n, r.instance) for r in rows))
    pairs = [("fit:triangle:mse", "fit:normal:mse"),
             ("fit:triangle:mse", "paper:empP"),
             ("fit:normal:mse", "paper:empP")]
    for (a, b) in pairs
        agree_mse = 0; agree_ew = 0; tot = 0
        for key in instances
            ra = get(bykey, (key..., a), nothing)
            rb = get(bykey, (key..., b), nothing)
            (ra === nothing || rb === nothing) && continue
            dr = ra.regret_nrm - rb.regret_nrm
            abs(dr) < 1e-12 && continue
            tot += 1
            sign(ra.mse - rb.mse) == sign(dr) && (agree_mse += 1)
            sign(ra.ew_avg - rb.ew_avg) == sign(dr) && (agree_ew += 1)
        end
        println(rpad("$a vs $b", 40),
                " mse: $(agree_mse)/$tot   ew_avg: $(agree_ew)/$tot")
    end

    # ---- Part C: weighted-norm fits vs entrywise fits ---------------------
    println("\n== Part C: fit in weighted norm vs entrywise (paired, regret_nrm) ==")
    for model in ["triangle", "normal"]
        wins = 0; ties = 0; losses = 0; diffs = Float64[]
        for key in instances
            rm = get(bykey, (key..., "fit:$model:mse"), nothing)
            rw = get(bykey, (key..., "fit:$model:wgt"), nothing)
            (rm === nothing || rw === nothing) && continue
            d = rm.regret_nrm - rw.regret_nrm   # >0 ⇒ weighted fit better
            push!(diffs, d)
            d > 1e-12 ? (wins += 1) : d < -1e-12 ? (losses += 1) : (ties += 1)
        end
        println(rpad(model, 10), "wgt better/tie/worse: $wins/$ties/$losses,  ",
                "mean Δregret (mse−wgt) = ", fmt(mean(diffs)))
    end

    # ---- raw vs normalized objective, by variant class --------------------
    println("\n== Normalize-don't-veto: regret_raw − regret_nrm by class ==")
    classes = [("exactN", r -> r.variant == "exactN"),
               ("perturbations", r -> startswith(r.variant, "pert:")),
               ("fits (all)", r -> startswith(r.variant, "fit:")),
               ("paper:empP", r -> r.variant == "paper:empP"),
               ("paper:binP", r -> r.variant == "paper:binP")]
    println(rpad("class", 16), "mean Δ     nrm better/tie/worse")
    for (name, pred) in classes
        sel = filter(pred, rows)
        d = [r.regret_raw - r.regret_nrm for r in sel]
        w = count(>(1e-12), d); t = count(x -> abs(x) <= 1e-12, d); l = count(<(-1e-12), d)
        println(rpad(name, 16), fmt(mean(d)), "    $w/$t/$l")
    end
end

main()
