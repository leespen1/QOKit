#=
E020 analysis — leakage vs K (per angle) and binned-proxy regret vs K.
Pure stdlib. Run: julia research/experiments/020_binned-weighted/analyze.jl
=#

using Statistics

function main()
    leak = Dict{Tuple{String, Int, String, Int}, Vector{Float64}}()   # (fam,n,angle,K)
    reg = Dict{Tuple{String, Int, Int, Int}, Vector{Float64}}()       # (fam,n,p,K)
    for f in sort(filter(f -> occursin(r"^results_task\d+\.csv$", f), readdir(@__DIR__)))
        for (i, line) in enumerate(eachline(joinpath(@__DIR__, f)))
            i == 1 && continue
            p = split(line, ",")
            fam, n, kind, K = String(p[1]), parse(Int, p[2]), String(p[6]), parse(Int, p[7])
            if startswith(kind, "leak")
                push!(get!(leak, (fam, n, kind, K), Float64[]), parse(Float64, p[9]))
            elseif kind == "binned_proxy"
                push!(get!(reg, (fam, n, parse(Int, p[8]), K), Float64[]),
                      parse(Float64, p[11]))
            end
        end
    end
    Ks = sort(unique(k[4] for k in keys(leak)))
    println("== one-layer leakage vs K (means) ==")
    println(rpad("cell / angle", 34), join(rpad.(string.(Ks), 8)))
    for key in sort(unique((k[1], k[2], k[3]) for k in keys(leak)))
        vals = [haskey(leak, (key..., K)) ?
                rpad(string(round(mean(leak[(key..., K)]); digits=4)), 8) :
                rpad("--", 8) for K in Ks]
        println(rpad("$(key[1]) n=$(key[2]) $(key[3])", 34), join(vals))
    end
    println("\n== binned-proxy regret vs K (means) ==")
    println(rpad("cell", 24), join(rpad.(string.(Ks), 8)))
    for key in sort(unique((k[1], k[2], k[3]) for k in keys(reg)))
        vals = [haskey(reg, (key..., K)) ?
                rpad(string(round(mean(reg[(key..., K)]); digits=4)), 8) :
                rpad("--", 8) for K in Ks]
        println(rpad("$(key[1]) n=$(key[2]) p=$(key[3])", 24), join(vals))
    end
end

main()
