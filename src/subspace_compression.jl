#=
subspace_compression.jl — The homogeneous proxy viewed as a subspace compression.

The homogeneous proxy tracks one amplitude per cost value. Equivalently, it
restricts the QAOA evolution to the subspace of statevectors that are constant
on each cost class S_c = {x : c(x) = c}. The orthogonal projection P onto this
subspace replaces each amplitude by its cost-class mean. One proxy step with
the same-instance empirical N(c'; d, c) is exactly: apply one true QAOA layer,
then project (in the per-bitstring "common amplitude" coordinates Q(c)).

The functions here implement that projection and the exactly-compressed
trajectory, plus the per-layer leakage ‖(I−P) U φ‖ that quantifies how much
of the state escapes the subspace at each layer. All operations are O(2^n)
per layer; no O(4^n) distribution and no 2^n × 2^n matrices are ever formed.
=#

export project_onto_cost_classes, reconstruct_from_cost_classes
export compressed_qaoa_trajectory


"""
    project_onto_cost_classes(state, costs; num_costs=maximum(costs)+1)

Orthogonally project a length-2^n statevector onto the cost-class subspace,
returning the compressed coordinates rather than the full projected vector.

Returns a NamedTuple:
- `Q`: length-`num_costs` vector; `Q[1+c]` is the mean amplitude over
  bitstrings of cost `c` (the homogeneous-proxy coordinates). Zero for
  unattained costs.
- `counts`: class sizes `M_c` (`counts[1+c] = |{x : c(x) = c}|`).
- `residual_norm`: `‖(I−P) state‖`, the part of `state` outside the subspace.

The projected vector itself is `reconstruct_from_cost_classes(Q, costs)`, and
its norm satisfies `‖P state‖² = Σ_c counts[c]·|Q[c]|² = ‖state‖² − residual_norm²`.
"""
function project_onto_cost_classes(
    state::AbstractVector,
    costs::AbstractVector{<:Real};
    num_costs::Integer=Int(maximum(costs)) + 1,
)
    @assert length(state) == length(costs) "state and costs must have the same length (2^n)"
    T = eltype(state)
    sums = zeros(T, num_costs)
    counts = zeros(Int, num_costs)
    @inbounds for x in eachindex(state, costs)
        ci = Int(costs[x]) + 1
        sums[ci] += state[x]
        counts[ci] += 1
    end
    Q = similar(sums)
    @inbounds for ci in 1:num_costs
        Q[ci] = counts[ci] > 0 ? sums[ci] / counts[ci] : zero(T)
    end
    residual2 = zero(real(T))
    @inbounds for x in eachindex(state, costs)
        residual2 += abs2(state[x] - Q[Int(costs[x]) + 1])
    end
    # Guard against tiny negative values from cancellation
    return (; Q, counts, residual_norm=sqrt(max(residual2, zero(residual2))))
end


"""
    reconstruct_from_cost_classes(Q, costs)

Expand compressed proxy coordinates back to a full statevector: the amplitude
at bitstring `x` is `Q[1 + c(x)]`. Inverse of the compression on the subspace;
composing with `project_onto_cost_classes` recovers `Q` exactly (for attained
costs).
"""
function reconstruct_from_cost_classes(Q::AbstractVector, costs::AbstractVector{<:Real})
    return [Q[Int(c) + 1] for c in costs]
end


"""
    compressed_qaoa_trajectory(costs, n, γs, βs; num_costs=maximum(costs)+1)

Evolve the exactly-compressed (homogeneous-proxy) QAOA trajectory by
alternating one true QAOA layer with orthogonal projection back onto the
cost-class subspace, alongside the uncompressed trajectory for comparison.
Mathematically identical to `QAOA_proxy_basic` run with the same-instance
empirical homogeneous distribution, but computed in O(p·n·2^n) without any
distribution array.

Angles are in raw radians (same convention as `qaoa_statevector`).

Returns a NamedTuple of per-layer diagnostics (index ℓ+1 ↔ layer ℓ, with
layer 0 = initial uniform state where applicable):
- `Qs`: compressed states `Q_0 … Q_p` (each a length-`num_costs` vector of
  per-bitstring class amplitudes).
- `counts`: cost-class sizes `M_c` (layer-independent).
- `leakage`: `λ_ℓ = ‖(I−P) U_ℓ φ_{ℓ-1}‖` for ℓ = 1…p, the norm lost to
  outside the subspace at each layer.
- `compressed_norm`: `‖φ_ℓ‖` for ℓ = 0…p; decays exactly by the accumulated
  leakage (`‖φ_ℓ‖² = ‖φ_{ℓ-1}‖² − λ_ℓ²`).
- `overlap`: `|⟨ψ_ℓ|φ_ℓ⟩|` for ℓ = 0…p against the true (unprojected)
  trajectory `ψ_ℓ`; note `φ_ℓ` is unnormalized.
- `distance`: `‖ψ_ℓ − φ_ℓ‖` for ℓ = 0…p (Theorem-2 bound target:
  `distance[1+p] ≤ sum(leakage)`).
"""
function compressed_qaoa_trajectory(
    costs::AbstractVector{T}, n::Integer,
    γs::AbstractVector{<:Real}, βs::AbstractVector{<:Real};
    num_costs::Integer=Int(maximum(costs)) + 1,
) where T <: AbstractFloat
    @assert length(γs) == length(βs) "γs and βs must have the same length"
    @assert length(costs) == 1 << n "costs must have length 2^n"
    p = length(γs)
    CT = Complex{T}

    # Uniform superposition lies inside the subspace, so both trajectories
    # start at the same state.
    ψ = fill(CT(1 / sqrt(T(1 << n))), 1 << n)   # true trajectory
    φ = copy(ψ)                                  # compressed trajectory

    proj0 = project_onto_cost_classes(φ, costs; num_costs)
    counts = proj0.counts

    compressed_norm2(Q) = sum(counts[ci] * abs2(Q[ci]) for ci in eachindex(Q))

    Qs = [proj0.Q]
    leakage = zeros(real(T), p)
    compressed_norm = zeros(real(T), p + 1)
    overlap = zeros(real(T), p + 1)
    distance = zeros(real(T), p + 1)

    compressed_norm[1] = sqrt(compressed_norm2(proj0.Q))
    overlap[1] = abs(sum(conj(a) * b for (a, b) in zip(ψ, φ)))
    distance[1] = sqrt(sum(abs2(a - b) for (a, b) in zip(ψ, φ)))

    # Costs as the statevector code expects them
    costsT = convert(AbstractVector{T}, costs)::AbstractVector{T}

    for ℓ in 1:p
        apply_phase_gate!(ψ, costsT, γs[ℓ])
        apply_x_mixer!(ψ, βs[ℓ], n)

        apply_phase_gate!(φ, costsT, γs[ℓ])
        apply_x_mixer!(φ, βs[ℓ], n)
        proj = project_onto_cost_classes(φ, costs; num_costs)
        leakage[ℓ] = proj.residual_norm
        push!(Qs, proj.Q)
        φ = reconstruct_from_cost_classes(proj.Q, costs)

        compressed_norm[ℓ + 1] = sqrt(compressed_norm2(proj.Q))
        overlap[ℓ + 1] = abs(sum(conj(a) * b for (a, b) in zip(ψ, φ)))
        distance[ℓ + 1] = sqrt(sum(abs2(a - b) for (a, b) in zip(ψ, φ)))
    end

    return (; Qs, counts, leakage, compressed_norm, overlap, distance)
end
