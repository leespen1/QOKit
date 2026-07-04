#=
E018 analysis — per-cell comparison of transfer regret (this experiment)
against the exact-compression and sampled-N proxy regrets recorded by
exps 002 (n=12,14) and 010 (n=16,18). Pure stdlib.

Run: julia research/experiments/018_transfer-calibration/analyze.jl
=#

using Statistics

function cellmean!(d, key, val)
    push!(get!(d, key, Float64[]), val)
end

function main()
    # transfer regrets, keyed (family, n, p); mean over sources and instances,
    # plus best-of-3-sources per instance
    transfer = Dict{Tuple{String, Int, Int}, Vector{Float64}}()
    bytarget = Dict{Tuple{String, Int, Int, Int}, Vector{Float64}}()
    for (i, line) in enumerate(eachline(joinpath(@__DIR__, "results.csv")))
        i == 1 && continue
        p = split(line, ",")
        fam, n, inst = String(p[1]), parse(Int, p[2]), parse(Int, p[3])
        depth, reg = parse(Int, p[6]), parse(Float64, p[9])
        cellmean!(transfer, (fam, n, depth), reg)
        cellmean!(bytarget, (fam, n, depth, inst), reg)
    end
    best = Dict{Tuple{String, Int, Int}, Vector{Float64}}()
    for ((fam, n, depth, _), regs) in bytarget
        cellmean!(best, (fam, n, depth), minimum(regs))
    end

    # proxy regrets from 002 and 010
    exact = Dict{Tuple{String, Int, Int}, Vector{Float64}}()
    samp = Dict{Tuple{String, Int, Int}, Vector{Float64}}()
    f002 = joinpath(@__DIR__, "..", "002_baselines-and-headroom", "results.csv")
    for (i, line) in enumerate(eachline(f002))
        i == 1 && continue
        p = split(line, ",")
        fam, n = String(p[1]), parse(Int, p[2])
        cellmean!(exact, (fam, n, 1), parse(Float64, p[9]) - parse(Float64, p[10]))
        cellmean!(exact, (fam, n, 3), parse(Float64, p[13]) - parse(Float64, p[14]))
    end
    dir010 = joinpath(@__DIR__, "..", "010_scaleup-ranking")
    for f in sort(filter(f -> occursin(r"^results_task\d+\.csv$", f), readdir(dir010)))
        for (i, line) in enumerate(eachline(joinpath(dir010, f)))
            i == 1 && continue
            p = split(line, ",")
            fam, n = String(p[1]), parse(Int, p[2])
            ex1 = parse(Float64, p[10]); ex3 = parse(Float64, p[13])
            isnan(ex1) || cellmean!(exact, (fam, n, 1), parse(Float64, p[8]) - ex1)
            isnan(ex3) || cellmean!(exact, (fam, n, 3), parse(Float64, p[11]) - ex3)
            cellmean!(samp, (fam, n, 1), parse(Float64, p[8]) - parse(Float64, p[9]))
            cellmean!(samp, (fam, n, 3), parse(Float64, p[11]) - parse(Float64, p[12]))
        end
    end

    fams = ["ER(0.5)", "ER(0.25)", "BA(k=2)", "BA(k=4)",
            "WS(k=4;b=0.1)", "WS(k=4;b=0.5)", "3-regular"]
    fmt(d, k) = haskey(d, k) ? rpad(string(round(mean(d[k]); digits=4)), 8) : rpad("--", 8)
    for depth in (1, 3)
        println("== p=$depth: mean regret per cell ==")
        println(rpad("family", 15), "n   transfer  best-of-3  exact-N   sampled-N")
        for fam in fams, n in (12, 14, 16, 18)
            k = (fam, n, depth)
            haskey(transfer, k) || continue
            println(rpad(fam, 15), rpad(n, 4),
                    fmt(transfer, k), "  ", fmt(best, k), "   ", fmt(exact, k), "  ", fmt(samp, k))
        end
        println()
    end

    t_all = vcat([transfer[k] for k in keys(transfer)]...)
    e_all = vcat([exact[k] for k in keys(exact)]...)
    println("pooled transfer regret: mean $(round(mean(t_all); digits=4)), ",
            "median $(round(median(t_all); digits=4)), max $(round(maximum(t_all); digits=4))")
    println("pooled exact-N proxy regret: mean $(round(mean(e_all); digits=4))")
end

main()
