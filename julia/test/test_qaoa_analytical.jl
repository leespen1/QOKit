#=
Analytical reference tests for QAOA simulation.

These tests verify the CPU (and optionally GPU) QAOA implementations against
closed-form solutions derived by hand for small cases, providing a ground-truth
check that does not depend on any external code (Python/QOKit).

═══════════════════════════════════════════════════════════════════════════════
Derivation: 2-qubit single-edge QAOA
═══════════════════════════════════════════════════════════════════════════════

Graph: n=2 vertices, single edge (0,1).
Costs: c(00)=0, c(01)=1, c(10)=1, c(11)=0  →  costs = [0, 1, 1, 0].

By symmetry, the QAOA state always has the form [u, v, v, u].

Initial state (uniform superposition): u₀ = v₀ = 1/2.

One QAOA layer (γ, β) evolves [u, v, v, u] as follows:

  Phase gate:  e^{-iγC/2} multiplies cost-1 amplitudes by φ = cis(-γ/2).
               State becomes [u, vφ, vφ, u].

  X-mixer:     e^{-iβX₀} ⊗ e^{-iβX₁} with c = cos(β), s = sin(β).
               After applying X₀ then X₁ to [u, vφ, vφ, u], the result is:
                 u_new = cos(2β)·u  - i·sin(2β)·v·φ
                 v_new = cos(2β)·v·φ - i·sin(2β)·u

  (Proof sketch: X₀ pairs indices (0,1) and (2,3). X₁ pairs (0,2) and (1,3).
   After X₀: [(cu - isvφ), (cvφ - isu), (cvφ - isu), (cu - isvφ)].
   After X₁ on pairs (0,2): new₀ = c(cu-isvφ) - is(cvφ-isu)
     = c²u - icsvφ - icsvφ - s²u = cos(2β)u - i·sin(2β)·vφ.
   Pair (1,3) gives the same by the [·,·,·,·] symmetry.)

For p=1, the expectation value has a closed form:
  ⟨C⟩ = 2|v|² = (1 + sin(4β)·sin(γ/2)) / 2
=#

using JuliaQAOA, Test


# ─── maxcut_costs ────────────────────────────────────────────────────────────

@testset "maxcut_costs on known graphs" begin
    # Single edge (0,1): cost = 1 iff bits differ
    @test maxcut_costs(2, [(0,1)]) == [0.0, 1.0, 1.0, 0.0]

    # Path P₃: edges (0,1), (1,2)
    # 000→0, 001→1, 010→2, 011→1, 100→1, 101→2, 110→1, 111→0
    @test maxcut_costs(3, [(0,1), (1,2)]) == [0, 1, 2, 1, 1, 2, 1, 0]

    # Triangle K₃: every non-trivial partition cuts exactly 2 of 3 edges
    @test maxcut_costs(3, [(0,1), (1,2), (0,2)]) == [0, 2, 2, 2, 2, 2, 2, 0]

    # No edges: all costs zero
    @test maxcut_costs(2, Tuple{Int,Int}[]) == [0, 0, 0, 0]
end


# ─── Phase gate ──────────────────────────────────────────────────────────────

@testset "phase gate on known state" begin
    # e^{-iγC/2}: each amplitude multiplied by cis(-γ·c(x)/2)
    costs = [0.0, 1.0, 1.0, 0.0]
    γ = 1.0
    state = ComplexF64[0.5, 0.3+0.1im, 0.2-0.4im, 0.1+0.2im]
    original = copy(state)
    apply_phase_gate!(state, costs, γ)

    φ = cis(-γ / 2)
    expected = [original[1], original[2]*φ, original[3]*φ, original[4]]
    @test state ≈ expected atol=1e-15
end


# ─── X-mixer ─────────────────────────────────────────────────────────────────

@testset "X-mixer from |00⟩" begin
    # e^{-iβX₀} ⊗ e^{-iβX₁} applied to |00⟩ gives:
    #   [cos²β, -i cosβ sinβ, -i cosβ sinβ, -sin²β]
    β = 0.7
    c, s = cos(β), sin(β)
    state = ComplexF64[1, 0, 0, 0]
    apply_x_mixer!(state, β, 2)

    expected = ComplexF64[c^2, -im*c*s, -im*c*s, -s^2]
    @test state ≈ expected atol=1e-15
end

@testset "X-mixer preserves norm" begin
    state = ComplexF64[0.3+0.1im, 0.5-0.2im, -0.1+0.4im, 0.2-0.3im]
    state ./= sqrt(sum(abs2, state))
    apply_x_mixer!(state, 1.234, 2)
    @test sum(abs2, state) ≈ 1.0 atol=1e-14
end


# ─── Full QAOA: 2-qubit single-edge analytical recurrence ───────────────────

"""
    analytical_2qubit_qaoa(γs, βs)

Compute the exact QAOA state for n=2, edge (0,1), using the recurrence:
  u_new = cos(2β)·u  - i·sin(2β)·v·cis(-γ/2)
  v_new = cos(2β)·v·cis(-γ/2) - i·sin(2β)·u

Returns (u, v) such that the full state is [u, v, v, u].
"""
function analytical_2qubit_qaoa(γs, βs)
    u = 0.5 + 0.0im
    v = 0.5 + 0.0im
    for ℓ in eachindex(γs)
        φ = cis(-γs[ℓ] / 2)
        c2β = cos(2βs[ℓ])
        s2β = sin(2βs[ℓ])
        u_new = c2β * u  - im * s2β * v * φ
        v_new = c2β * v * φ - im * s2β * u
        u, v = u_new, v_new
    end
    return u, v
end

@testset "QAOA p=1 matches analytical" begin
    γ, β = 1.23, 0.456
    costs = maxcut_costs(2, [(0,1)])

    u, v = analytical_2qubit_qaoa([γ], [β])

    state = qaoa_statevector(costs, 2, [γ], [β])
    @test state ≈ [u, v, v, u] atol=1e-14

    exp_val = qaoa_expectation(costs, 2, [γ], [β])
    @test exp_val ≈ 2 * abs2(v) atol=1e-14
end

@testset "QAOA p=3 matches analytical recurrence" begin
    γs = [0.8, 1.5, 0.3]
    βs = [0.4, 0.9, 0.2]
    costs = maxcut_costs(2, [(0,1)])

    u, v = analytical_2qubit_qaoa(γs, βs)

    state = qaoa_statevector(costs, 2, γs, βs)
    @test state ≈ [u, v, v, u] atol=1e-13

    exp_val = qaoa_expectation(costs, 2, γs, βs)
    @test exp_val ≈ 2 * abs2(v) atol=1e-13
end


# ─── p=1 expectation closed form ─────────────────────────────────────────────

@testset "p=1 expectation closed form: ⟨C⟩ = (1 + sin(4β)sin(γ/2))/2" begin
    costs = maxcut_costs(2, [(0,1)])
    for (γ, β) in [(1.0, 0.5), (π/3, π/4), (2.0, 0.1), (0.01, 1.5)]
        exp_val = qaoa_expectation(costs, 2, [γ], [β])
        analytical = (1 + sin(4β) * sin(γ/2)) / 2
        @test exp_val ≈ analytical atol=1e-14
    end
end
