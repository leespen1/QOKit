################################################################################
# 
# Trying out the trasnfer matrix approach to doing QAOA. I want to default to
# this when given a homogeneous distribution as an array (instead of being
# given the proxy).
#
# Later on, I can adapt it to GPU.
#
################################################################################

"""
Build the transfer matrix (in-place) which, when multiplied onto the vector 
[exp(-im*gamma*c)*Q(c) for c in 0:num_constraints], will perform a QAOA proxy
layer update (i.e. Equation (5) in the paper).

IMPORTANT! Ordering is N(c', c, d), for maximum efficiency.

For now, it is assumed that the ordering of the homogeneous distribution is
N(c',d,c). In the future, I think it would be more efficient to order it as
N(c',c,d).
"""
function build_transfer_matrix!(
    M::AbstractArray{ComplexF64, 2}, 
    homodist::AbstractArray{<: Real, 3},
    beta::Real,
)
    @assert size(homodist,1) == size(homodist,2) "Dimensions 1 and 2 of homodist must be the same length."
    @assert size(M,1) == size(M,2) == size(homodist,1) "Transfer matrix must be square, and size consistent with homogeneous distribution."

    num_costs = size(homodist, 1)
    num_distances = size(homodist, 3)
    num_qubits = num_distances - 1


    d_vals = 0:num_qubits
    cosb = cos(beta)
    neg1_im_sinb = -1im*sin(beta)
    weights = @. (cosb ^ (num_qubits - d_vals)) * (neg1_im_sinb ^ d_vals)
    M .= 0
    for i in 1:num_distances
        homodist_slice = view(homodist, :, :, i)
        @. M += weights[i] * homodist_slice
    end
    return M
end

"""
Allocating version of `build_transfer_matrix!`
"""
function build_transfer_matrix(
    homodist::AbstractArray{<: Real, 3},
    beta::Real,
)
    num_costs = size(homodist, 1)
    M = Matrix{ComplexF64}(undef, num_costs, num_costs)
    build_transfer_matrix!(M, homodist, beta)
    return M
end



"""
Apply one layer of the QAOA proxy circuit
- `state_vec::AbstractVector{ComplexF64}`: The proxy state vector at the end of
the previous layer.
- `transfer_matrix::AbstractMatrix{ComplexF64}`: The transfer matrix for this
layer.
- `gamma::Real`: The value of gamma at this layer.
- `const_vals::AbstractVector{Int64}`: [0, ..., num_costs-1]  (maybe could
replace with 0:num_costs-1, not sure if it will allocate memory)
- `state_vec_storage::AbstractVector{ComplexF64}`: A vector with the same
length as state_vec.
"""
function qaoa_proxy_circuit_layer!(
    state_vec::AbstractVector{<: Complex},
    transfer_matrix::AbstractMatrix{<: Complex},
    gamma::Real,
    cost_vals::AbstractVector{<: Integer},
    state_vec_storage::AbstractVector{<: Complex},
)
    @. state_vec_storage = exp(-1im * gamma * cost_vals)
    state_vec_storage .*= state_vec
    mul!(state_vec, transfer_matrix, state_vec_storage)
    return state_vec 
end

function qaoa_proxy_circuit(
    homodist::AbstractArray{<: Real, 3},
    gammas::AbstractVector{<: Real},
    betas::AbstractVector{<: Real},
)
    @assert size(homodist,1) == size(homodist,2) "Dimensions 1 and 2 of homodist must be the same length."
    @assert length(gammas) == length(betas) "Gamma and beta vectors must be the same length."

    num_layers = length(gammas)
    num_costs = size(homodist, 1)
    cost_vals = collect(0:num_costs-1)
    num_distances = size(homodist, 3)
    num_qubits = num_distances - 1

    state_mat = zeros(ComplexF64, num_costs, 1+num_layers)
    state_vec = zeros(ComplexF64, num_costs)
    state_vec_storage = zeros(ComplexF64, num_costs)

    init_amplitude = sqrt(1 / (1 << num_qubits))
    state_vec .= init_amplitude
    state_mat[:,1] .= state_vec


    transfer_matrix = zeros(ComplexF64, num_costs, num_costs)
    for i in 1:num_layers
        beta = betas[i]
        gamma = gammas[i]
        build_transfer_matrix!(transfer_matrix, homodist, beta)
        qaoa_proxy_circuit_layer!(state_vec, transfer_matrix, gamma, cost_vals, state_vec_storage)

        state_mat[:,1+i] .= state_vec
    end

    return state_mat
end

function qaoa_proxy_circuit(
    proxy,
    gammas::AbstractVector{<: Real},
    betas::AbstractVector{<: Real},
)
    N = efficient_order_homodist(proxy)
    return qaoa_proxy_circuit(N, gammas, betas)
end

