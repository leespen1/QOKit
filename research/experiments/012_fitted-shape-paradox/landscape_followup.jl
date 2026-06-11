#=
E2.4 follow-up: amplitude error proved to be the wrong norm too (see
analysis_summary.txt — regret tracks neither mse_entry nor ‖Q₁_model−Q₁_emp‖).
Regret is determined by *where the model landscape puts its argmax*, so the
candidate predictors are landscape similarity and argmax displacement, both
invariant to the scale conventions that wreck the amplitude norm.

Rebuilds every model N from the params recorded in results.csv (no refitting),
re-evaluates the p=1 proxy landscape on the same 40×40 grid, and records per
(instance, model, raw/cal):
  pearson, spearman — correlation of E_model(γ,β) with E_emp(γ,β) over the grid
  dgamma, dbeta     — argmax displacement from the empirical proxy's argmax
  same_argmax       — whether the grid argmax coincides with the empirical one
  slicesum_min/max  — range of N(c';:,:) slice sums over attained classes
                      (documents each model's scale convention; emp = 2^n)
Writes landscape.csv and prints Spearman correlations of regret (joined from
results.csv) against the landscape metrics.

Run from the repo root (~5 min with threads):
  JULIA_NUM_THREADS=auto julia --project research/experiments/012_fitted-shape-paradox/landscape_followup.jl
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads
using Statistics: mean, median, std

const DIR = @__DIR__
const SEED = 20260611
const NS = [12, 14]
const INSTANCES = 10
const P1_GRID_LEN = 40
const P1_γ = collect(range(0.0, π; length=P1_GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=P1_GRID_LEN))

const FAMILIES = [
    ("ER(0.5)",       (rng, n) -> erdos_renyi_edges(n, 0.5; rng)),
    ("ER(0.25)",      (rng, n) -> erdos_renyi_edges(n, 0.25; rng)),
    ("BA(k=2)",       (rng, n) -> barabasi_albert_edges(n, 2; rng)),
    ("BA(k=4)",       (rng, n) -> barabasi_albert_edges(n, 4; rng)),
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("WS(k=4;b=0.5)", (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]

instance_seed(fam_idx, n, inst) = SEED + 10_000 * fam_idx + 100 * n + inst

function calibrate_slices(N::Array{Float64,3}, n::Integer)
    out = copy(N)
    for i in axes(out, 1)
        s = sum(@view out[i, :, :])
        s > 0 && (@view(out[i, :, :]) .*= (2.0^n / s))
    end
    return out
end

function tiedrank(v)
    p = sortperm(v)
    r = zeros(length(v))
    i = 1
    while i <= length(v)
        j = i
        while j < length(v) && v[p[j+1]] == v[p[i]]
            j += 1
        end
        r[p[i:j]] .= (i + j) / 2
        i = j + 1
    end
    return r
end

function pearson(x, y)
    sx, sy = std(x), std(y)
    (sx == 0 || sy == 0) && return NaN
    return mean((x .- mean(x)) .* (y .- mean(y))) / (sx * sy)
end

spearman(x, y) = pearson(tiedrank(x), tiedrank(y))

function main()
    # results.csv supplies the fitted params and the regrets to correlate against
    lines = readlines(joinpath(DIR, "results.csv"))
    header = split(lines[1], ',')
    col = Dict(h => i for (i, h) in enumerate(header))
    recs = [split(l, ',') for l in lines[2:end] if !isempty(strip(l))]
    params = Dict((r[col["family"]], parse(Int, r[col["n"]]),
                   parse(Int, r[col["instance"]]), r[col["model"]]) => r[col["params"]]
                  for r in recs)
    regret = Dict((r[col["family"]], parse(Int, r[col["n"]]),
                   parse(Int, r[col["instance"]]), r[col["model"]], tag) =>
                      parse(Float64, r[col["regret_" * tag]])
                  for r in recs, tag in ("raw", "cal"))

    jobs = []
    for (fam_idx, (fam_name, gen)) in enumerate(FAMILIES), n in NS, inst in 1:INSTANCES
        seed = instance_seed(fam_idx, n, inst)
        edges = gen(MersenneTwister(seed), n)
        push!(jobs, (; fam_name, n, inst, edges))
    end

    homodists = Vector{Array{Float64,3}}(undef, length(jobs))
    for k in eachindex(jobs)
        j = jobs[k]
        costs = maxcut_costs(j.n, j.edges)
        homodists[k] = get_homogeneous_distribution_from_costs_direct(
            costs, length(j.edges), j.n)
    end
    println("empirical homodists built")

    grid_pts = vec([(γ, β) for γ in P1_γ, β in P1_β])
    γmat = reshape([g for (g, _) in grid_pts], :, 1)
    βmat = reshape([b for (_, b) in grid_pts], :, 1)

    getp(s) = parse.(Float64, split(s, ';'))

    results = Vector{String}(undef, length(jobs))
    @threads for k in eachindex(jobs)
        j = jobs[k]
        n, m = j.n, length(j.edges)
        costs = maxcut_costs(n, j.edges)
        N_emp = homodists[k]
        counts = zeros(Int, m + 1)
        for cst in costs
            counts[Int(cst) + 1] += 1
        end
        P_emp = counts ./ (1 << n)

        Q_emp = QAOA_proxy_multi(N_emp, γmat, βmat)[end]
        E_emp = vec(expectation(Q_emp, P_emp, n))
        i_emp = argmax(E_emp)

        key(mdl) = (j.fam_name, n, j.inst, mdl)
        models = [
            ("emp", N_emp),
            ("paper", cpu_compute_homodist(PaperProxy(m, n, 2m / (n * (n - 1))))),
            [(mdl, cpu_compute_homodist(OldTriangleProxy(m, n, getp(params[key(mdl)])...)))
             for mdl in ("tri_init", "tri_25", "tri_50", "tri_fit")]...,
            [(mdl, cpu_compute_homodist(NormalProxy(m, n, getp(params[key(mdl)])...)))
             for mdl in ("norm_init", "norm_25", "norm_50", "norm_fit")]...,
        ]

        rows = String[]
        for (name, N_model) in models
            sums = [sum(@view N_model[i, :, :]) for i in axes(N_model, 1)]
            attained = [s for (i, s) in enumerate(sums) if sum(@view N_emp[i, :, :]) > 0]
            for (tag, N_used) in (("raw", N_model), ("cal", calibrate_slices(N_model, n)))
                Q = QAOA_proxy_multi(N_used, γmat, βmat)[end]
                E = vec(expectation(Q, P_emp, n))
                i_mod = argmax(E)
                push!(rows, join(Any[j.fam_name, n, j.inst, name, tag,
                                     pearson(E, E_emp), spearman(E, E_emp),
                                     abs(grid_pts[i_mod][1] - grid_pts[i_emp][1]),
                                     abs(grid_pts[i_mod][2] - grid_pts[i_emp][2]),
                                     Int(i_mod == i_emp),
                                     minimum(attained), maximum(attained)], ","))
            end
        end
        results[k] = join(rows, "\n")
        println("done: $(j.fam_name) n=$n inst=$(j.inst)")
    end

    outpath = joinpath(DIR, "landscape.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,model,tag,pearson,spearman_grid," *
                    "dgamma,dbeta,same_argmax,slicesum_min,slicesum_max")
        foreach(r -> println(io, r), results)
    end
    println("landscape metrics → $outpath\n")

    # --- correlate regret with the landscape metrics -------------------------
    lrows = NamedTuple[]
    for l in readlines(outpath)[2:end]
        f = split(l, ',')
        push!(lrows, (family = f[1], n = parse(Int, f[2]), instance = parse(Int, f[3]),
                      model = String(f[4]), tag = String(f[5]),
                      pearson = parse(Float64, f[6]), sp = parse(Float64, f[7]),
                      dγ = parse(Float64, f[8]), dβ = parse(Float64, f[9]),
                      same = parse(Int, f[10])))
    end
    MODELS = ["paper", "tri_init", "tri_25", "tri_50", "tri_fit",
              "norm_init", "norm_25", "norm_50", "norm_fit"]
    instances = unique((r.family, r.n, r.instance) for r in lrows)
    byim = Dict((r.family, r.n, r.instance, r.model, r.tag) => r for r in lrows)
    for tag in ("raw", "cal")
        for (label, xf) in [("1 − pearson(landscape)", r -> 1 - r.pearson),
                            ("argmax displacement ‖(Δγ,Δβ)‖", r -> hypot(r.dγ, r.dβ))]
            ρs = Float64[]
            for kk in instances
                xs = [xf(byim[(kk..., mdl, tag)]) for mdl in MODELS]
                ys = [regret[(kk..., mdl, tag)] for mdl in MODELS]
                ρ = spearman(xs, ys)
                isnan(ρ) || push!(ρs, ρ)
            end
            println("regret_$tag ~ $label:  per-instance Spearman " *
                    "$(round(mean(ρs), digits=3)) ± $(round(std(ρs), digits=3)), " *
                    "med $(round(median(ρs), digits=3))  (N=$(length(ρs)))")
        end
        same_frac = Dict(mdl => mean(byim[(kk..., mdl, tag)].same for kk in instances)
                         for mdl in ["emp"; MODELS])
        println("  argmax coincides with emp ($tag): ",
                join(("$mdl=$(round(same_frac[mdl], digits=2))" for mdl in ["emp"; MODELS]), " "))
    end
end

main()
