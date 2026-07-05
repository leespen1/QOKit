#=
E027 — Does the map survive a different problem? Max-3-XOR replication.

Theorems 1–3 hold for any cost function; every measurement so far is
MaxCut. Max-3-XOR (each clause: x_i ⊕ x_j ⊕ x_k = b, cost = #satisfied)
is sud2024's second family, with integer costs — native cost classes, no
binning. Replicated here at two clause densities:

  (i)  one-layer leakage at three reference angles (density-law analogue);
  (ii) parameter setting at p=1 (40×40, γ ∈ [0, 2π]) and p=3 ramps
       (γ ∈ [0.05, 3.2]): exactN, sampledN (S=10), transfer×3 from n=10
       same-density sources, pooled universal — against true grid ceilings.

Cells: m ∈ {4n, 8n} × n ∈ {14, 16}, 10 instances (tasks 1–4).
Submit:  cd research/experiments/027_max3xor && sbatch run.sb
Smoke:   E27_SMOKE=1 julia --project research/experiments/027_max3xor/run.jl
=#

using JuliaQAOA
using Random: MersenneTwister, rand, randperm
using Base.Threads: @threads
using Statistics: mean

const SMOKE = get(ENV, "E27_SMOKE", "0") == "1"

const SEED = 20260718
const SRC_SEED = 20260719
const SOURCE_N = 10
const N_SOURCES = SMOKE ? 1 : 3
const S_N = 10
const NINST = SMOKE ? 2 : 10
const GRID_LEN = SMOKE ? 10 : 40
const RAMP_LEN = SMOKE ? 4 : 8
const P_RAMP = 3
const LEAK_ANGLES = [(0.2, 0.2), (0.5, 0.3), (1.0, 0.4)]

const P1_γ = collect(range(0.0, 2π; length=GRID_LEN))
const P1_β = collect(range(0.0, π/2; length=GRID_LEN))
const RAMP_γ = collect(range(0.05, 3.2; length=RAMP_LEN))
const RAMP_β = collect(range(0.05, 0.8; length=RAMP_LEN))

const TASKS = [("4n", 14), ("8n", 14), ("4n", 16), ("8n", 16)]
mclauses(dens, n) = dens == "4n" ? 4n : 8n

"Random 3-XOR instance: (i,j,k,b) clauses, distinct vars per clause."
function random_3xor(rng, n, m)
    clauses = Vector{NTuple{4, Int}}(undef, m)
    for c in 1:m
        v = randperm(rng, n)[1:3]
        clauses[c] = (v[1] - 1, v[2] - 1, v[3] - 1, rand(rng, 0:1))
    end
    clauses
end

"Costs: number of satisfied parity clauses, threaded."
function xor_costs(n, clauses)
    costs = zeros(Float64, 1 << n)
    @threads for x in 0:((1 << n) - 1)
        c = 0
        for (i, j, k, b) in clauses
            c += (((x >> i) & 1) ⊻ ((x >> j) & 1) ⊻ ((x >> k) & 1)) == b
        end
        costs[x + 1] = c
    end
    costs
end

"One-layer leakage from |+⟩ onto the native cost classes."
function one_layer_leakage(costs, m, n, γ, β)
    ψ = fill(ComplexF64(1 / sqrt(2.0^n)), 1 << n)
    apply_phase_gate!(ψ, costs, γ)
    apply_x_mixer!(ψ, β, n)
    sums = zeros(ComplexF64, m + 1)
    sizes = zeros(Int, m + 1)
    for i in eachindex(ψ)
        b = Int(costs[i]) + 1
        sums[b] += ψ[i]
        sizes[b] += 1
    end
    means = [sizes[b] > 0 ? sums[b] / sizes[b] : zero(ComplexF64) for b in 1:(m + 1)]
    resid2 = 0.0
    for i in eachindex(ψ)
        resid2 += abs2(ψ[i] - means[Int(costs[i]) + 1])
    end
    sqrt(resid2)
end

function true_vals(costs, n, scheds)
    vals = zeros(length(scheds))
    @threads for k in eachindex(scheds)
        γs, βs = scheds[k]
        vals[k] = qaoa_expectation(costs, n, γs, βs)
    end
    vals
end

function main()
    tid = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
    dens, n = SMOKE ? ("4n", 12) : TASKS[tid == 0 ? 1 : tid]
    println("task $tid: density=$dens n=$n")

    p1_scheds = vec([([γ], [β]) for γ in P1_γ, β in P1_β])
    ramp_combos = vec([(g1, gf, b1, bf) for g1 in RAMP_γ, gf in RAMP_γ,
                                            b1 in RAMP_β, bf in RAMP_β])
    ramp_scheds = [linear_ramp(c..., P_RAMP) for c in ramp_combos]

    # sources at n=10, same clause density
    srcv = Vector{Dict{Int, Vector{Float64}}}()
    for rep in 1:N_SOURCES
        srng = MersenneTwister(SRC_SEED + rep + 100 * (dens == "8n"))
        sclauses = random_3xor(srng, SOURCE_N, mclauses(dens, SOURCE_N))
        scosts = xor_costs(SOURCE_N, sclauses)
        c_opt = maximum(scosts)
        push!(srcv, Dict(1 => true_vals(scosts, SOURCE_N, p1_scheds) ./ c_opt,
                         3 => true_vals(scosts, SOURCE_N, ramp_scheds) ./ c_opt))
    end
    transfer_idx = Dict(p => [argmax(srcv[r][p]) for r in 1:N_SOURCES] for p in (1, 3))
    universal_idx = Dict(p => argmax(reduce(+, [r[p] for r in srcv])) for p in (1, 3))

    rows = String[]
    for inst in 1:NINST
        seed = SEED + 100 * n + 1000 * (dens == "8n") + inst
        rng = MersenneTwister(seed)
        m = mclauses(dens, n)
        clauses = random_3xor(rng, n, m)
        costs = xor_costs(n, clauses)
        c_opt = maximum(costs)

        for (γ, β) in LEAK_ANGLES
            λ = one_layer_leakage(costs, m, n, γ, β)
            push!(rows, join(Any[dens, n, inst, seed, m, "leak_g$(γ)_b$(β)", 1,
                                 λ, "", ""], ","))
        end

        counts = zeros(Int, m + 1)
        for c in costs
            counts[Int(c) + 1] += 1
        end
        P_emp = counts ./ (1 << n)
        N_exact = get_homogeneous_distribution_from_costs_direct(costs, m, n)
        N_samp = sampled_homogeneous_distribution(costs, m, n;
                     samples_per_class=S_N, rng=MersenneTwister(seed + 777))

        for (p, scheds, γmat_βmat) in ((1, p1_scheds, nothing), (P_RAMP, ramp_scheds, nothing))
            tv = true_vals(costs, n, scheds)
            ceil_ar = maximum(tv) / c_opt
            emit(method, idx) = begin
                ar = tv[idx] / c_opt
                @assert ceil_ar - ar > -1e-8
                push!(rows, join(Any[dens, n, inst, seed, m, method, p, ar,
                                     ceil_ar, ceil_ar - ar], ","))
            end
            push!(rows, join(Any[dens, n, inst, seed, m, "ceiling", p, ceil_ar,
                                 ceil_ar, 0.0], ","))
            if p == 1
                γmat = reshape([s[1][1] for s in scheds], :, 1)
                βmat = reshape([s[2][1] for s in scheds], :, 1)
            else
                γmat, βmat = linear_ramp_matrix(
                    [c[1] for c in ramp_combos], [c[2] for c in ramp_combos],
                    [c[3] for c in ramp_combos], [c[4] for c in ramp_combos], P_RAMP)
            end
            for (name, N) in (("exactN", N_exact), ("sampledN", N_samp))
                Qs = QAOA_proxy_multi(N, γmat, βmat)
                emit(name, argmax(vec(expectation(Qs[end], P_emp, n))))
            end
            for r in 1:N_SOURCES
                emit("transfer_$r", transfer_idx[p][r])
            end
            emit("universal", universal_idx[p])
        end
        println("done: inst=$inst m=$m")
        flush(stdout)
    end

    suffix = SMOKE ? "_smoke" : "_task$tid"
    outpath = joinpath(@__DIR__, "results$suffix.csv")
    open(outpath, "w") do io
        println(io, "density,n,instance,seed,m,method,p,ar,ceil,regret")
        foreach(r -> println(io, r), rows)
    end
    println("E027 task complete → $outpath")
end

main()
