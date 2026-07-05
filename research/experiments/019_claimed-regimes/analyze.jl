#=
E019 analysis — per-cell mean regret by method. Pure stdlib.
Run: julia research/experiments/019_claimed-regimes/analyze.jl
=#

using Statistics

function main()
    cells = Dict{Tuple{String, Int, Int}, Dict{String, Vector{Float64}}}()
    for f in sort(filter(f -> occursin(r"^results_task\d+\.csv$", f), readdir(@__DIR__)))
        for (i, line) in enumerate(eachline(joinpath(@__DIR__, f)))
            i == 1 && continue
            p = split(line, ",")
            key = (String(p[1]), parse(Int, p[2]), parse(Int, p[6]))
            d = get!(cells, key, Dict{String, Vector{Float64}}())
            push!(get!(d, String(p[7]), Float64[]), parse(Float64, p[10]))
        end
    end
    methods = ["exactN", "sampledN", "paper_binP", "paper_empP",
               "transfer_1", "transfer_2", "transfer_3", "universal"]
    hdr = ["exactN", "sampled", "pap_binP", "pap_empP", "tr1", "tr2", "tr3", "univ"]
    println(rpad("cell", 24), join(rpad.(hdr, 9)))
    for key in sort(collect(keys(cells)); by=k -> (k[3] <= 3 ? 1 : 0, k[2], k[1], k[3]))
        fam, n, p = key
        d = cells[key]
        vals = [haskey(d, m) ? rpad(string(round(mean(d[m]); digits=4)), 9) :
                rpad("--", 9) for m in methods]
        println(rpad("$fam n=$n p=$p", 24), join(vals))
    end
    println("\nbest transfer (min over the 3 sources, per instance) per cell:")
    for key in sort(collect(keys(cells)); by=k -> (k[3] <= 3 ? 1 : 0, k[2], k[1], k[3]))
        d = cells[key]
        all(haskey(d, "transfer_$r") for r in 1:3) || continue
        ni = length(d["transfer_1"])
        best = [minimum(d["transfer_$r"][i] for r in 1:3) for i in 1:ni]
        println(rpad("$(key[1]) n=$(key[2]) p=$(key[3])", 24),
                round(mean(best); digits=4))
    end
end

main()
