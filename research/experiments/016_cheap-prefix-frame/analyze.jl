#=
Analysis for E016. Reads results.csv and, per ramp regime, compares three
frames for capturing the full depth-P=20 trajectory:
  E_cc      cost-class frame, dim m+1, ZERO statevector layers (the proxy's frame)
  E_full_d  optimal top-d PCA of the full trajectory (oracle ceiling, needs all P)
  E_prefix_d  top-d PCA of only the cheap K-layer prefix (the constructive method)
and the largest principal angle between the prefix-d and full-d subspaces
(subspace-discoverability metric; 0° = the prefix already spans the late subspace).

Run: julia --project research/experiments/016_cheap-prefix-frame/analyze.jl
=#
using Statistics: mean

const DIR = @__DIR__
function readcsv(path)
    lines = readlines(path); h = String.(split(lines[1], ','))
    cols = Dict(c => String[] for c in h)
    for ln in lines[2:end]; isempty(ln) && continue; f = split(ln, ',')
        for (j, c) in enumerate(h); push!(cols[c], String(f[j])); end; end
    cols
end
fc(c, k) = parse.(Float64, c[k])

c = readcsv(joinpath(DIR, "results.csv"))
ramp = c["ramp"]; Ecc = fc(c,"E_cc"); Eful = fc(c,"E_full_d"); Epre = fc(c,"E_prefix_d")
ang = fc(c,"max_angle_deg")
d = Int(fc(c,"d")[1]); K = Int(fc(c,"K")[1]); P = Int(fc(c,"P")[1])

println("Frames for the full P=$P trajectory:  cost-class (dim m+1, 0 layers) vs")
println("optimal top-$d PCA (oracle, needs all $P) vs cheap top-$d PCA of the first $K layers.\n")
println(rpad("ramp",10), rpad("E_cc",9), rpad("E_full_$d",10), rpad("E_pre_$d",9),
        rpad("max∠(pre,full)",16), "prefix>cost-class?")
for r in ["small","moderate","large"]
    s = ramp .== r
    frac = count(Epre[s] .> Ecc[s]) / count(s)
    println(rpad(r,10), rpad(round(mean(Ecc[s]),digits=4),9),
            rpad(round(mean(Eful[s]),digits=4),10), rpad(round(mean(Epre[s]),digits=4),9),
            rpad(string(round(mean(ang[s]),digits=1),"°"),16),
            "$(round(100*frac))% of instances")
end
println("\nReadings:")
println("• E_full_$d ≈ $(round(mean(Eful),digits=3)) confirms the trajectory is ~$d-dim (exp 009).")
println("• Cheap-prefix capture E_pre vs cost-class E_cc, and the principal angle, say")
println("  whether the late subspace is cheaply discoverable from an early prefix.")
