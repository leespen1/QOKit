#=
E016 — The Discussion's "most interesting open question": the QAOA ramp
trajectory lives in a 2–4-dim moving subspace (exp 009), but the proxy uses the
fixed (m+1)-dim cost-class frame and pays for the mis-aim. Could an
instance-adapted low-rank frame beat it — and can such a frame be built CHEAPLY,
breaking the chicken-and-egg of "needing states to build the frame"?

Concrete test. For a depth-P ramp, build a d-dim frame from only a cheap PREFIX
(the first K+1 statevector layers, K≪P), then ask:
  - subspace stability: principal angles between the prefix's top-d PCA subspace
    and the FULL trajectory's top-d PCA subspace. Small ⇒ the early trajectory
    already spans the late subspace ⇒ the frame is cheaply discoverable.
  - capture: fraction of the FULL trajectory's energy captured by the prefix
    frame (E_prefix_d), vs the optimal full-PCA frame (E_full_d, oracle ceiling)
    and vs the cost-class frame (E_cc, dim m+1, the proxy's frame — built from
    costs alone, zero statevector layers).

Honesty: the prefix frame costs K statevector layers (O(K·2ⁿ)) — cheaper than the
full P but NOT free, and exponential like simulation. This is a STRUCTURAL test
(is the right subspace cheaply discoverable / stable over depth), not a speedup
claim. The cost-class frame remains the only zero-layer option.

Outputs results.csv: family,n,inst,seed,m,ramp,P,K,d,E_cc,E_full_d,E_prefix_d,
max_angle_deg (largest principal angle, prefix-d vs full-d).

Run (CPU, threaded; ~5 min):
  JULIA_NUM_THREADS=auto julia --project research/experiments/016_cheap-prefix-frame/run.jl
Smoke test: E16_SMOKE=1 julia --project .../run.jl
=#

using JuliaQAOA
using Random: MersenneTwister
using Base.Threads: @threads
using LinearAlgebra: svd, svdvals

const SMOKE = get(ENV, "E16_SMOKE", "0") == "1"
const SEED  = 20260611
const NS    = SMOKE ? [10] : [12, 14]
const INST  = SMOKE ? 2 : 5
const P     = 20
const D     = 4          # target frame dimension (≈ k99 from exp 009)
const KPRE  = 5          # prefix length: layers 0..KPRE (KPRE+1 states), ≈ P/4

const RAMPS = [
    (id="small",    g1=0.1, gf=0.4, b1=0.4, bf=0.1),
    (id="moderate", g1=0.1, gf=0.8, b1=0.6, bf=0.1),
    (id="large",    g1=0.2, gf=1.6, b1=0.8, bf=0.1),
]
const FAMILIES = [
    ("ER(0.5)",       (rng, n) -> erdos_renyi_edges(n, 0.5; rng)),
    ("ER(0.25)",      (rng, n) -> erdos_renyi_edges(n, 0.25; rng)),
    ("BA(k=2)",       (rng, n) -> barabasi_albert_edges(n, 2; rng)),
    ("BA(k=4)",       (rng, n) -> barabasi_albert_edges(n, 4; rng)),
    ("WS(k=4;b=0.1)", (rng, n) -> watts_strogatz_edges(n, 4, 0.1; rng)),
    ("WS(k=4;b=0.5)", (rng, n) -> watts_strogatz_edges(n, 4, 0.5; rng)),
    ("3-regular",     (rng, n) -> random_regular_edges(n, 3; rng)),
]

# fraction of trajectory energy captured by the orthonormal-column frame U
function captured(U, states)
    tot = 0.0
    for ψ in states
        c = U' * ψ                  # coordinates in the frame
        tot += sum(abs2, c)         # ‖U Uᴴ ψ‖² = ‖Uᴴ ψ‖²  (U has orthonormal cols)
    end
    return tot / length(states)     # since ‖ψ‖²=1
end

function main()
    jobs = [(fi, fam, gen, n, inst, ramp)
            for (fi, (fam, gen)) in enumerate(FAMILIES)
            for n in NS for inst in 1:INST for ramp in RAMPS]
    res = Vector{String}(undef, length(jobs))
    @threads for k in eachindex(jobs)
        fi, fam, gen, n, inst, ramp = jobs[k]
        seed = SEED + 10_000 * fi + 100 * n + inst
        edges = gen(MersenneTwister(seed), n); m = length(edges)
        costs = maxcut_costs(n, edges)
        γs, βs = linear_ramp(ramp.g1, ramp.gf, ramp.b1, ramp.bf, P)
        states = qaoa_statevector(costs, n, γs, βs; return_intermediates = true)

        # cost-class frame (dim m+1): captured energy from the projection residual
        Ecc = 0.0
        for ψ in states
            r = project_onto_cost_classes(ψ, costs; num_costs = m + 1).residual_norm
            Ecc += 1 - r^2
        end
        Ecc /= length(states)

        # full-trajectory top-D PCA frame (oracle ceiling)
        Ufull = svd(stack(states)).U[:, 1:D]
        Efull = captured(Ufull, states)

        # cheap-prefix top-D PCA frame
        Upre = svd(stack(states[1:KPRE+1])).U[:, 1:D]
        Epre = captured(Upre, states)

        # principal angles between prefix-D and full-D subspaces
        cosang = svdvals(Upre' * Ufull); clamp!(cosang, 0.0, 1.0)
        max_angle = acosd(minimum(cosang))

        res[k] = join((fam, n, inst, seed, m, ramp.id, P, KPRE, D,
            round(Ecc, digits=5), round(Efull, digits=5), round(Epre, digits=5),
            round(max_angle, digits=2)), ",")
        println("done: ", rpad(fam,14), " n=$n inst=$inst ", rpad(ramp.id,9),
                " Ecc=$(round(Ecc,digits=3)) Efull=$(round(Efull,digits=3)) ",
                "Epre=$(round(Epre,digits=3)) angle=$(round(max_angle,digits=1))°")
    end
    path = joinpath(@__DIR__, SMOKE ? "results_smoke.csv" : "results.csv")
    open(path, "w") do io
        println(io, "family,n,inst,seed,m,ramp,P,K,d,E_cc,E_full_d,E_prefix_d,max_angle_deg")
        foreach(r -> println(io, r), res)
    end
    println("\nwrote $path  ($(length(jobs)) jobs)")
end

main()
