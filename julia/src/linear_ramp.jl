"""
    linear_ramp(γ₁, γ_f, β₁, β_f, p; pi_units=true)

Generate a linear ramp QAOA parameter schedule with p layers.

The paper (Eq. 13-14) defines:
    γ_ℓ = γ₁ + (γ_f - γ₁) · ℓ/p    for ℓ = 1, …, p
    β_ℓ = β₁ + (β_f - β₁) · ℓ/p    for ℓ = 1, …, p

This reduces the 2p-dimensional parameter space to just 4 parameters,
making high-depth (large p) optimization tractable.

Returns (γs, βs) as vectors of length p.

If pi_units=true (default), parameters are in units of π (consistent with
QAOA_proxy_single/multi). If false, parameters are in raw radians.
"""
function linear_ramp(γ₁::Real, γ_f::Real, β₁::Real, β_f::Real, p::Int; pi_units::Bool=true)
    layers = collect(1:p) ./ p
    γs = γ₁ .+ (γ_f - γ₁) .* layers
    βs = β₁ .+ (β_f - β₁) .* layers
    return γs, βs
end

"""
    linear_ramp_matrix(γ₁s, γ_fs, β₁s, β_fs, p; pi_units=true)

Generate linear ramp schedules for multiple parameter sets at once.
Each input is a vector of length K (number of parameter sets).

Returns (γs, βs) as K×p matrices, compatible with QAOA_proxy_multi.
"""
function linear_ramp_matrix(
    γ₁s::AbstractVector{<:Real}, γ_fs::AbstractVector{<:Real},
    β₁s::AbstractVector{<:Real}, β_fs::AbstractVector{<:Real},
    p::Int; pi_units::Bool=true
)
    K = length(γ₁s)
    @assert length(γ_fs) == K && length(β₁s) == K && length(β_fs) == K
    layers = collect(1:p)' ./ p  # 1×p row vector
    γs = γ₁s .+ (γ_fs .- γ₁s) .* layers  # K×p
    βs = β₁s .+ (β_fs .- β₁s) .* layers  # K×p
    return γs, βs
end
