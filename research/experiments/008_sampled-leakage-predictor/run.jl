#=
E2.2 — Can per-class leakage be estimated from a few sampled bitstrings per
cost class?

Theorem 3: the one-layer leakage of class state |c'⟩ satisfies
    λ(c')² = M_{c'}^{-1} Σ_c M_c · Var_{y∈S_c}[ g_{c'}(y) ],
with g_{c'}(y) = Σ_d f_d(β) n(y; d, c'), f_d(β) = cos(β)^{n-d}(-i·sinβ)^d.
(The e^{-iγc'/2} phase is unimodular and drops out of the variance.)

If the within-class variances can be estimated from S sampled bitstrings per
class, leakage — the paper's central diagnostic — becomes computable at scales
where the exact O(4^n) distribution is out of reach. This experiment measures
the estimator's accuracy vs S, re-verifying old hypothesis H-sample in the
quantity that matters (it was originally observed for entrywise N estimation).

Note Var here is the *population* variance over the class; the estimator uses
the unbiased sample variance from S draws (without replacement when S ≥ M_c).

For each instance and a few (γ, β) points:
  - exact λ(c') for every attained class via one statevector layer + projection;
  - sampled λ̂(c') for S ∈ {2, 5, 10, 25} via stratified sampling;
  - report relative errors of λ̂ and of the aggregate η_F = sqrt(Σ λ(c')²).

This also constitutes the first machine-precision-level check of Theorem 3's
identity itself: at S = "all" the estimator must reproduce the exact λ(c').

Run from the repo root (~30 min with threads):
  JULIA_NUM_THREADS=auto julia --project research/experiments/008_sampled-leakage-predictor/run.jl
Smoke test: E1_SMOKE=1 julia --project ...
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads
using Statistics: mean, var

const SMOKE = get(ENV, "E1_SMOKE", "0") == "1"

const SEED = 20260611
const NS = SMOKE ? [10] : [12, 14]
const INSTANCES = SMOKE ? 2 : 10
const SS = [2, 5, 10, 25]
const ANGLES = [(0.2, 0.15), (0.4, 0.3), (0.8, 0.6)]   # small → moderate (γ, β)

const FAMILIES = [
    ("ER(0.5)",       (rng, n) -> erdos_renyi_edges(n, 0.5; rng)),
    ("ER(0.25)",      (rng, n) -> erdos_renyi_edges(n, 0.25; rng)),
    ("BA(k=2)",       (rng, n) -> barabasi_albert_edges(n, 2; rng)),
    ("BA(k=4)",       (rng, n) -> barabasi_albert_edges(n, 4; rng)),
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("WS(k=4;b=0.5)", (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]

"n(y; d, c) for one bitstring y: one O(2^n) scan."
function neighborhood_histogram(y::Int, costs, n, num_costs)
    h = zeros(Float64, n + 1, num_costs)
    @inbounds for z in 0:(length(costs) - 1)
        h[count_ones(y ⊻ z) + 1, Int(costs[z + 1]) + 1] += 1.0
    end
    return h
end

function main()
    jobs = [(fam_idx, fam_name, gen, n, inst)
            for (fam_idx, (fam_name, gen)) in enumerate(FAMILIES)
            for n in NS, inst in 1:INSTANCES]
    results = Vector{String}(undef, length(jobs))

    @threads for k in eachindex(jobs)
        fam_idx, fam_name, gen, n, inst = jobs[k]
        seed = SEED + 10_000 * fam_idx + 100 * n + inst
        rng = MersenneTwister(seed)
        edges = gen(rng, n)
        m = length(edges)
        costs = maxcut_costs(n, edges)
        num_costs = m + 1

        counts = zeros(Int, num_costs)
        members = [Int[] for _ in 1:num_costs]
        for x in 0:(length(costs) - 1)
            ci = Int(costs[x + 1]) + 1
            counts[ci] += 1
            push!(members[ci], x)
        end
        attained = findall(>(0), counts)

        work = zeros(ComplexF64, length(costs))
        rows = String[]
        for (γ, β) in ANGLES
            f_d = [cos(β)^(n - d) * (-im * sin(β))^d for d in 0:n]

            # Exact λ(c') via statevector layer + projection
            λ_exact = Dict{Int,Float64}()
            for ci in attained
                fill!(work, 0)
                for x in members[ci]
                    work[x + 1] = 1 / sqrt(counts[ci])
                end
                apply_phase_gate!(work, costs, γ)
                apply_x_mixer!(work, β, n)
                λ_exact[ci] = project_onto_cost_classes(work, costs; num_costs).residual_norm
            end
            η_exact = sqrt(sum(v^2 for v in values(λ_exact)))

            # Sampled estimator: g_{c'}(y) for sampled y, per class.
            # S = 0 means full enumeration (population variance) — at that
            # setting λ̂ must reproduce λ_exact, i.e. it verifies Theorem 3.
            S_list = n <= 12 ? vcat(SS, 0) : SS
            for S in S_list
                λ2̂ = zeros(num_costs)   # accumulates Σ_c M_c·Var, indexed by c'
                for ci in attained
                    Mc = counts[ci]
                    nsamp = S == 0 ? Mc : min(S, Mc)
                    ys = nsamp == Mc ? members[ci] :
                         members[ci][randperm_partial(rng, Mc, nsamp)]
                    gs = Matrix{ComplexF64}(undef, nsamp, num_costs)
                    for (si, y) in enumerate(ys)
                        h = neighborhood_histogram(y, costs, n, num_costs)
                        for cp in 1:num_costs
                            gs[si, cp] = sum(f_d[d + 1] * h[d + 1, cp] for d in 0:n)
                        end
                    end
                    if nsamp == 1
                        continue   # singleton sample: population variance is 0
                    end
                    # Population variance when the class is fully enumerated,
                    # unbiased sample variance otherwise
                    for cp in attained
                        v = var(view(gs, :, cp); corrected=(nsamp < Mc))
                        λ2̂[cp] += Mc * real(v)
                    end
                end
                for ci in attained
                    λ2̂[ci] /= counts[ci]
                end
                η̂ = sqrt(sum(λ2̂[cp] for cp in attained))
                relerr_η = abs(η̂ - η_exact) / η_exact
                med_rel = let errs = [abs(sqrt(max(λ2̂[ci], 0)) - λ_exact[ci]) /
                                      max(λ_exact[ci], 1e-12) for ci in attained]
                    sort(errs)[cld(length(errs), 2)]
                end
                push!(rows, join([fam_name, n, inst, seed, m, γ, β, S,
                                  η_exact, η̂, relerr_η, med_rel], ","))
            end
        end
        results[k] = join(rows, "\n")
        println("done: $fam_name n=$n inst=$inst")
    end

    outpath = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(outpath, "w") do io
        println(io, "family,n,instance,seed,m,gamma,beta,S," *
                    "eta_exact,eta_sampled,relerr_eta,median_relerr_lambda")
        foreach(r -> println(io, r), results)
    end
    println("\nE2.2 sweep complete → $outpath")
end

"Indices of `k` distinct uniform draws from 1:N (partial Fisher-Yates)."
function randperm_partial(rng, N, k)
    idx = collect(1:N)
    for i in 1:k
        j = rand(rng, i:N)
        idx[i], idx[j] = idx[j], idx[i]
    end
    return idx[1:k]
end

main()
