#=
Diagnostic for the PaperProxy anomaly found in experiment 004: on dense ER(0.5)
the analytical-N proxy's grid argmax lands on an unphysical peak (predicted
⟨C⟩ = 93 on a 38-edge graph) where the proxy state's norm has inflated ~7.5×.

Exact compression (Theorem 1) is contractive — the proxy norm can only decay —
so norm inflation is a zero-cost certificate of model error. Filtering grid
points to ‖φ‖² ≤ 1.05 restores PaperProxy to competitive parameters on the
same instance. Experiment 005 applies this filter across the full instance set.

Run from the repo root:
  julia --project research/experiments/004_gate-leakage-vs-regret/diagnose_paper_artifact.jl
=#

using JuliaQAOA
using Random: MersenneTwister

seed = 20260611 + 10_000 * 1 + 100 * 12 + 1   # ER(0.5), n=12, instance 1 (004 scheme)
n = 12
edges = erdos_renyi_edges(n, 0.5; rng=MersenneTwister(seed))
m = length(edges)
costs = maxcut_costs(n, edges)
c_opt = maximum(costs)

N_emp = get_homogeneous_distribution_from_costs_direct(costs, m, n)
counts = zeros(Int, m + 1)
for c in costs
    counts[Int(c) + 1] += 1
end
P_emp = counts ./ 2^n

pp = PaperProxy(m, n, 2m / (n * (n - 1)))
N_pap = cpu_compute_homodist(pp)
P_pap = [P_cost_distribution(pp, c) for c in 0:m]

combos = vec([(γ, β) for γ in range(0.0, π; length=40), β in range(0.0, π/2; length=40)])
γmat = reshape([c[1] for c in combos], :, 1)
βmat = reshape([c[2] for c in combos], :, 1)

proxy_vals(N, P) = vec(expectation(QAOA_proxy_multi(N, γmat, βmat)[end], P, n))
proxy_norm2s(N, P) = vec(2.0^n .* sum(abs2.(QAOA_proxy_multi(N, γmat, βmat)[end]) .* P, dims=1))
real_at(γβ) = qaoa_expectation(costs, n, [γβ[1]], [γβ[2]])

vals_emp = proxy_vals(N_emp, P_emp)
vals_pap = proxy_vals(N_pap, P_pap)
norm2_pap = proxy_norm2s(N_pap, P_pap)

i_e, i_p = argmax(vals_emp), argmax(vals_pap)
println("instance: ER(0.5) n=$n m=$m c_opt=$c_opt")
println("empirical-N argmax (γ,β)=", round.(combos[i_e], digits=3),
        ": real ⟨C⟩ = ", round(real_at(combos[i_e]), digits=2))
println("PaperProxy raw argmax (γ,β)=", round.(combos[i_p], digits=3),
        ": predicted ⟨C⟩ = ", round(vals_pap[i_p], digits=2),
        " (max possible: $(m)!), norm² = ", round(norm2_pap[i_p], digits=2),
        ", real ⟨C⟩ = ", round(real_at(combos[i_p]), digits=2))
@assert vals_pap[i_p] > m "expected the unphysical predicted ⟨C⟩ > m artifact"
@assert norm2_pap[i_p] > 2 "expected norm inflation at the artifact"

filtered = copy(vals_pap)
filtered[norm2_pap .> 1.05] .= -Inf
i_f = argmax(filtered)
println("PaperProxy norm-filtered argmax (γ,β)=", round.(combos[i_f], digits=3),
        ": real ⟨C⟩ = ", round(real_at(combos[i_f]), digits=2),
        "  [grid ceiling: ", round(maximum(real_at(c) for c in combos), digits=2), "]")
