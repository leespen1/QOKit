#=
figure4_squared_overlap.jl — Reproduce Figure 4 from the paper.

Paper Figure 4: Squared overlap |⟨ψ_true|ψ_proxy⟩|² vs QAOA layer ℓ,
for a single graph instance with linear ramp schedule at various (γ₁, γ_f).

The "proxy state" here is NOT the compressed proxy (which only tracks m+1
amplitudes). Instead, it keeps all 2^n amplitudes but substitutes
N(c(x); d, c) for n(x; d, c) in the evolution equation. This lets us
compute a meaningful overlap with the true statevector.

Paper parameters: G(8, 1/2), p=20, linear ramp, multiple (γ₁, γ_f) curves
Quick test:       G(6, 0.5), p=8, 3 curves

Started:  2026-03-17
Finished: 2026-03-17
=#

include("common.jl")

#==============================================================================#
#                          CONFIGURATION                                        #
#==============================================================================#

const N_QUBITS = 8
const P_EDGE = 0.5
const P_DEPTH = 20
const SEED = 42

# Linear ramp parameters: β₁ and β_f are shared across curves.
# Each curve uses a different (γ₁, γ_f) pair.
const β₁ = 0.45    # Starting β (in radians)
const β_f = 0.05   # Final β

# (γ₁, γ_f) pairs for different curves (paper uses several small-γ pairs)
# Note: γ uses QOKit convention where phase gate is exp(-iγC/2)
const GAMMA_PAIRS = [
    (0.10, 0.30),
    (0.20, 0.60),
    (0.30, 1.00),
    (0.50, 1.40),
]

# Which proxy to use for the homogeneous approximation
PROXY_CONFIG = ("PaperProxy", (m, n) -> PaperProxy(m, n, P_EDGE))

#==============================================================================#
#                    PROXY STATEVECTOR SIMULATION                               #
#==============================================================================#

"""
    proxy_statevector(costs, n, homodist, γs, βs)

Run the "extended proxy" that keeps all 2^n amplitudes but uses the
homogeneous distribution N(c(x); d, c) in place of n(x; d, c).

For each layer ℓ:
  q_ℓ(x) = Σ_{d,c} cos(β)^(n-d) (-i sin β)^d exp(-iγc) q_{ℓ-1}(y_c) N(c(x); d, c)

where the sum runs over all (d, c) pairs and y_c indexes bitstrings with
cost c (but since we use N instead of n, we weight by the homogeneous avg).

Actually, the correct interpretation is simpler: the proxy replaces
n(x; d, c) → N(c(x); d, c) in the evolution. So the update for bitstring
x at layer ℓ is:

  q_ℓ(x) = Σ_y [cos β]^(n-d(x,y)) [-i sin β]^d(x,y) exp(-iγ c(y)) q_{ℓ-1}(y)

with the approximation that N(c'; d, c) ≈ n(x; d, c) for all x with c(x)=c'.

Under this approximation, q_ℓ(x) depends on x only through c(x). So we:
1. Evolve the compressed proxy Q_ℓ(c') using the standard proxy algorithm
2. Reconstruct the full state: q_ℓ(x) = Q_ℓ(c(x)) for all x

Returns a vector of full statevectors (one per layer, including layer 0).
"""
function proxy_statevector_from_compressed(
    costs::Vector{Float64}, n::Int, m::Int,
    homodist::AbstractArray{<:Real, 3},
    γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real};
    pi_units::Bool=false
)
    p = length(γs)
    num_bitstrings = 1 << n

    # Convert to pi_units for the proxy if needed
    if pi_units
        γs_pi = γs
        βs_pi = βs
    else
        γs_pi = γs ./ π
        βs_pi = βs ./ π
    end

    # Run compressed proxy (returns list of Q vectors, length m+1 each)
    Qs = QAOA_proxy_single(homodist, γs_pi, βs_pi; pi_units=true)

    # Reconstruct full states from compressed Q
    states = Vector{Vector{ComplexF64}}(undef, p + 1)
    for ℓ in 0:p
        Q = Qs[ℓ + 1]
        state = Vector{ComplexF64}(undef, num_bitstrings)
        for x in 0:(num_bitstrings - 1)
            c = Int(costs[x + 1])
            # Clamp to valid range (costs may exceed proxy's m)
            c_clamped = min(c, m)
            state[x + 1] = Q[c_clamped + 1]
        end
        states[ℓ + 1] = state
    end

    return states
end

"""
    squared_overlap(state1, state2)

Compute |⟨state1|state2⟩|².
"""
function squared_overlap(state1::Vector{ComplexF64}, state2::Vector{ComplexF64})
    return abs2(dot(state1, state2))
end

#==============================================================================#
#                          COMPUTATION                                          #
#==============================================================================#

println("=== Figure 4: Squared Overlap vs QAOA Layer ===")
println("Parameters: n=$N_QUBITS, p_edge=$P_EDGE, p=$P_DEPTH")

# Generate one graph instance
rng = MersenneTwister(SEED)
inst = generate_er_instance(N_QUBITS, P_EDGE; rng)
m = inst.num_edges
println("Graph: $(inst.num_vertices) vertices, $m edges")

# Compute proxy distribution from the chosen proxy
label, constructor = PROXY_CONFIG
proxy = constructor(m, N_QUBITS)
homodist = cpu_compute_homodist(proxy)

# For each (γ₁, γ_f) pair, compute overlaps at each layer
println("Computing overlaps for $(length(GAMMA_PAIRS)) parameter curves...")
overlap_curves = map(GAMMA_PAIRS) do (γ1, γf)
    # Generate linear ramp parameters (in radians)
    γs, βs = linear_ramp(γ1, γf, β₁, β_f, P_DEPTH)

    # Real QAOA: evolve statevector layer by layer (GPU if available)
    real_states = qaoa_statevector_intermediates_device(inst.costs, N_QUBITS, γs, βs)

    # Proxy: evolve via compressed proxy, reconstruct full states
    proxy_states = proxy_statevector_from_compressed(
        inst.costs, N_QUBITS, m, homodist, γs, βs; pi_units=false
    )

    # Compute squared overlap at each layer
    overlaps = [squared_overlap(real_states[ℓ+1], proxy_states[ℓ+1]) for ℓ in 0:P_DEPTH]

    (γ1=γ1, γf=γf, overlaps=overlaps)
end

#==============================================================================#
#                          PLOTTING                                             #
#==============================================================================#

println("Plotting...")
fig = Figure(size=FIGURE_SIZE)
ax = Axis(fig[1, 1],
    xlabel="QAOA Layer ℓ",
    ylabel="Squared Overlap |⟨ψ_true|ψ_proxy⟩|²",
    title="Proxy vs True QAOA Overlap\nG($N_QUBITS, $P_EDGE), p=$P_DEPTH, $label",
)

colors = [:blue, :red, :green, :orange, :purple, :brown]
for (i, curve) in enumerate(overlap_curves)
    scatterlines!(ax, 0:P_DEPTH, curve.overlaps,
        label="γ₁=$(curve.γ1), γ_f=$(curve.γf)",
        color=colors[mod1(i, length(colors))],
        markersize=6,
    )
end

axislegend(ax, position=:lb)
ylims!(ax, 0, 1.05)

save_figure(fig, "figure4_squared_overlap.png")
println("Done!")
