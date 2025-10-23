#=
This file implements the interface for approximating QAOA using a homogenous
proxy. The main functions are:
- QAOA_proxy: *runs* the homogeneous proxy of the QAOA circuit.
- QAOA_proxy_expectation: computes the *expectation value* of a proxy state.

There is no equivalent to the python function `QAOA_proxy_optimize_gamma_beta`.
Instead, `QAOA_proxy_optimize_gamma_beta` should be called in python, but `proxy`
should be a Julia object in order to use these Julia functions throughout the
optimization.
=#

"""
New version, using homodist array, but not "fused matrix-multiplication."

Optionally, provide state_vec to avoid allocating the state_vec
"""
function QAOA_proxy(
    homodist::AbstractArray{<: Real, 3}, gammas::AbstractVector{<: Real},
    betas::AbstractVector{<: Real},
    state_vec1::AbstractVector{ComplexF64} = zeros(ComplexF64, size(homodist, 1)),
    state_vec2::AbstractVector{ComplexF64} = zeros(ComplexF64, size(homodist, 1)),
)::Vector{ComplexF64}
    @assert length(gammas) == length(betas) "Gamma vec and beta vec must be same length."
    @assert size(homodist, 1) == size(homodist, 3) "1st and 3rd dimensions of homogeneous distribution must be the same length."
    @assert length(state_vec1) == length(state_vec2) == size(homodist, 1) "State vector arrays must have some length as number of unique costs in homogeneous distribution."
    p = length(gammas)

    num_distances = size(homodist, 2)
    num_costs = size(homodist, 1)

    init_amplitude = 1 / sqrt(2 ^ (num_distances-1))
    state_vec1 .= init_amplitude

    costs = 0:num_costs-1
    distances = 0:num_distances-1
    β_factors = Vector{ComplexF64}(undef, num_distances)
    γ_factors = Vector{ComplexF64}(undef, num_costs)

    for ℓ in 1:p
        # If I do this with turbo, it's SIMD, so GPU translation should be easy (parallel across sets of β and γs is harder)
        sinβ, cosβ = sincos(betas[ℓ])
        neg_im_sinβ = -1im*sinβ
        γ = gammas[ℓ]

        map!(d -> cosβ^(num_distances-1-d) * neg_im_sinβ^(d), β_factors, distances)
        map!(c -> exp(-im*γ*c), γ_factors, costs)
        state_vec2 .= 0
        for i_c_prime in axes(homodist, 1), i_d in axes(homodist, 2), i_c in axes(homodist, 3)
            state_vec2[i_c_prime] += β_factors[i_d] * γ_factors[i_c] *
                                     state_vec1[i_c] * homodist[i_c_prime, i_d, i_c]
        end
        # swap state_vec1/2
        state_vec_tmp = state_vec1
        state_vec1 = state_vec2
        state_vec2 = state_vec_tmp
    end
    return state_vec1
end

"""
Version which uses a "fused matvec mult" approach. Should at least get some
speedup based on LienarAlgebra multithreading, but I'm not sure about the
actual algorithm. (Probably speedup will become more substantial when doing
multiple β's and γ's, due to superior time locality for mat-mat mults).
"""
function QAOA_proxy_matvec(
    homodist::AbstractArray{<: Real, 3}, gammas::AbstractVector{<: Real},
    betas::AbstractVector{<: Real},
    state_vec1::AbstractVector{ComplexF64} = zeros(ComplexF64, size(homodist, 1)),
    state_vec2::AbstractVector{ComplexF64} = zeros(ComplexF64, size(homodist, 1));
    use_BLAS=true,
)::Vector{ComplexF64}
    @assert length(gammas) == length(betas) "Gamma vec and beta vec must be same length."
    @assert size(homodist, 1) == size(homodist, 3) "1st and 3rd dimensions of homogeneous distribution must be the same length."
    @assert length(state_vec1) == length(state_vec2) == size(homodist, 1) "State vector arrays must have some length as number of unique costs in homogeneous distribution."
    p = length(gammas)

    num_distances = size(homodist, 2)
    num_costs = size(homodist, 1)

    init_amplitude = 1 / sqrt(2 ^ (num_distances-1))
    state_vec1 .= init_amplitude

    costs = 0:num_costs-1
    distances = 0:num_distances-1
    β_factors = Vector{ComplexF64}(undef, num_distances)
    γ_factors = Matrix{ComplexF64}(undef, 1, num_costs)
    # I think I can construct the vector v using broadcasted multiplication
    # with a column vector and row vector, then reshaping into a vector.

    M = reshape(homodist, size(homodist,1), :) # Reshape to wide matrix
    v_mat = Matrix{ComplexF64}(undef, num_distances, num_costs)
    v_vec_view = reshape(v_mat, :)


    for ℓ in 1:p
        # If I do this with turbo, it's SIMD, so GPU translation should be easy (parallel across sets of β and γs is harder)
        sinβ, cosβ = sincos(betas[ℓ])
        neg_im_sinβ = -1im*sinβ
        γ = gammas[ℓ]

        map!(d -> cosβ^(num_distances-1-d) * neg_im_sinβ^(d), β_factors, distances)
        map!(c -> exp(-im*γ*c), γ_factors, costs)

        state_row_vec = reshape(state_vec1, 1, :)
        @. v_mat = γ_factors * state_row_vec * β_factors

        if use_BLAS
            mul!(state_vec2, M, v_vec_view)
        else
            state_vec2 .= M*v_vec_view 
        end

        # swap state_vec1/2
        state_vec_tmp = state_vec1
        state_vec1 = state_vec2
        state_vec2 = state_vec_tmp
    end
    return state_vec1
end

"""
Version which uses a "fused matmat mult" approach, for multiple sets of gamma
and beta.
"""
function QAOA_proxy_matmat(
    homodist::AbstractArray{<: Real, 3}, gammas::AbstractVecOrMat{<: Real},
    betas::AbstractVecOrMat{<: Real},
    state_vecs1::AbstractMatrix{<: Complex} = zeros(ComplexF64, size(homodist, 1), size(gammas, 2)),
    state_vecs2::AbstractMatrix{<: Complex} = zeros(ComplexF64, size(homodist, 1), size(gammas, 2));
    use_BLAS=true
)
    @assert size(gammas) == size(betas) "Gamma vec and beta vec must be same shape."
    @assert size(homodist, 1) == size(homodist, 3) "1st and 3rd dimensions of homogeneous distribution must be the same length."
    @assert size(state_vecs1, 1) == size(state_vecs2, 1) == size(homodist, 1) "State vector arrays must have some length as number of unique costs in homogeneous distribution."
    @assert size(state_vecs1, 2) == size(state_vecs2, 2) == size(gammas, 2) == size(betas, 2) "Batch size must be consistent across gammas, betas, and statevecs."

    p = size(gammas, 1)
    num_batches = size(gammas, 2)
    num_distances = size(homodist, 2)
    num_costs = size(homodist, 1)

    init_amplitude = 1 / sqrt(2 ^ (num_distances-1))
    state_vecs1 .= init_amplitude

    costs = similar(homodist, Int32, 1, num_costs)
    copyto!(costs, 0:num_costs-1)
    distances = similar(homodist, Int32, num_distances)
    distances .= 0:num_distances-1

    β_factors = similar(state_vecs1, num_distances, 1, num_batches)
    γ_factors = similar(state_vecs1, 1, num_costs, num_batches)
    # I think I can construct the vector v using broadcasted multiplication
    # with a column vector and row vector, then reshaping into a vector.

    M = reshape(homodist, size(homodist,1), :) # Reshape to wide matrix
    v_3D = similar(state_vecs1, num_distances, num_costs, num_batches)
    v_mat_view = reshape(v_3D, :, num_batches)

    for ℓ in 1:p
        sinβ = reshape(sin.(betas[ℓ,:]), 1, 1, :)
        cosβ = reshape(cos.(betas[ℓ,:]), 1, 1, :)
        neg_im_sinβ = -1im .* sinβ
        γ = reshape(gammas[ℓ,:], 1, 1, :)

        @. β_factors = cosβ^(num_distances-1-distances) * neg_im_sinβ^(distances)
        @. γ_factors = exp(-im*γ*costs)

        state_row_vecs = reshape(state_vecs1, 1, :, num_batches)
        @. v_3D = γ_factors * state_row_vecs * β_factors

        if use_BLAS # For some reason, mul! seems to invoke generic matmul, instead of BLAS gemm
            # Problem seems to be that M is real, while the others are complex.
            # Interestingly, BLAS is smart enough to use gemv when M is real
            # and the others are vectors, but not if the others are matrices.
            # Also note that 3-arg mul uses 5-arg mul under the hood, and LinearAlgebra.BLAS.gemv! can be called directly if desired.
            @show typeof(state_vecs2), typeof(M), typeof(v_mat_view)
            #println(@code_typed mul!(state_vecs2, M, v_mat_view))
            mul!(state_vecs2, M, v_mat_view)
        else
            #println(@code_typed M*v_mat_view)
            state_vecs2 .= M*v_mat_view 
        end


        # swap state_vec1/2
        state_vecs_tmp = state_vecs1
        state_vecs1 = state_vecs2
        state_vecs2 = state_vecs_tmp
    end
    return state_vecs1
end


"""
Run the homogeneous proxy of the QAOA circuit with parameters gammas (cost) 
and betas (mixer).

Return a complex matrix where each row is the set of proxy amplitudes
after each "layer" of the QAOA circuit (i.e., the homogeneous proxy of the 
state vector).
"""
function QAOA_proxy(
        proxy,
        gammas::AbstractVector{<: Real},
        betas::AbstractVector{<: Real}
    )::Matrix{ComplexF64}

    @assert length(gammas) == length(betas)

    init_amplitude = sqrt(1 / (1 << proxy.num_qubits))

    num_QAOA_layers = length(gammas)
    num_costs = proxy.num_constraints + 1

    proxy_amplitudes = zeros(ComplexF64, num_costs, 1+num_QAOA_layers)
    proxy_amplitudes[:,1] .= init_amplitude
    
    for current_depth in 1:num_QAOA_layers
        prev_amplitudes = view(proxy_amplitudes, :, current_depth)
        gamma = gammas[current_depth]
        beta = betas[current_depth]

        for cost_1 in 0:num_costs-1
            proxy_amplitudes[1+cost_1, 1+current_depth] = compute_amplitude_sum(
                proxy, prev_amplitudes, gamma, beta, cost_1
            )
        end
    end

    return proxy_amplitudes
end

"""
Convert numpy arrays to julia arrays before doing QAOA_proxy.
(should check whether this improves performance, I don't think it should much)
"""
function QAOA_proxy(proxy, gamma::PyArray, beta::PyArray)::Matrix{ComplexF64}
    return QAOA_proxy(
        proxy,
        pyconvert(Vector, gamma),
        pyconvert(Vector, beta),
    )
end

"""
Given a set of proxy amplitudes (i.e. the homomgeneous proxy of the state
vector.), compute the expectation value associated with those amplitudes.
"""
function QAOA_proxy_expectation(
        proxy,
        proxy_amplitudes::AbstractVector{<: Number},
        drop_first_N_costs::Integer=0
    )::Float64

    proxy_expectation_value::Float64 = 0.0
    num_unique_costs = length(proxy_amplitudes)

    for cost in drop_first_N_costs:num_unique_costs-1
        proxy_expectation_value += cost * N_cost_distribution(proxy, cost) * abs2(proxy_amplitudes[1+cost])
    end

    return proxy_expectation_value
end

"""
Given parameters gamma and beta, and the set of probability amplitudes
Q_l(c) for each unique cost c after performing l "layers" of the QAOA
homogenous proxy, and a particular cost c', return Q_l+1(c'), which is the
amplitude of cost c' after the l+1-th layer of the QAOA homogeneous proxy.

See the for-loop of Algorithm 1 in parameter-setting paper.
"""
function compute_amplitude_sum(
        proxy,
        prev_amplitudes::AbstractVector{ComplexF64},
        gamma::Real,
        beta::Real,
        cost_1::Integer
    )::ComplexF64

    @assert length(prev_amplitudes) == proxy.num_constraints + 1

    gamma_factors = exp.(-1im * gamma * (0:proxy.num_constraints))
    d_vals = 0:proxy.num_qubits
    sinb, cosb = sincos(beta)
    neg1_im_sinb = -1im * sinb
    beta_factors = @. (cosb ^ (proxy.num_qubits - d_vals)) * (neg1_im_sinb ^ d_vals)
    sum::ComplexF64 = 0
    # Changed loop order for efficiency. Check for same result
    for distance in 0:proxy.num_qubits 
        beta_factor = beta_factors[1+distance]
        for cost_2 in 0:proxy.num_constraints
            gamma_factor = gamma_factors[1+cost_2]
            num_costs_at_distance = N_cost_distance_distribution(
                proxy, cost_1, distance, cost_2
            )
            sum += beta_factor * gamma_factor * prev_amplitudes[1+cost_2] * num_costs_at_distance
        end
    end
    return sum
end




"""
Convert numpy arrays to julia arrays before doing QAOA_proxy_expectation_python.
(should check whether this impacts performance)
"""
function QAOA_proxy_expectation(
        proxy,
        proxy_amplitudes::PyArray,
        drop_first_N_costs::Integer=0
    )::Float64
    @assert ndims(proxy_amplitudes) == 1
    return QAOA_proxy_expectation(
        proxy,
        pyconvert(Vector, proxy_amplitudes),
        drop_first_N_costs
    )
end

"""
Inverse here is the additive inverse, not multiplicative inverse.

I.e. we are taking -x, not 1/x.

Scipy allows us to minimize functions, but we want to maximize our cost
(the number of edges cut), so we have to flip the +/- sign to turn our
problem into a minimization problem.
"""
function inverse_proxy_objective_function(
        proxy,
        num_QAOA_layers::Integer,
        expectations::Union{Nothing, AbstractVector}=nothing
    )

    function inverse_objective(args...)::Float64
        # Note that args[0] in Python is args[1] in Julia 
        gammas_betas_combined_vec = pyconvert(Vector, args[1]) # Conversion may not be necessary. Remove if impacts performance significantly.

        @assert length(gammas_betas_combined_vec) == 2*num_QAOA_layers
        gammas = @view gammas_betas_combined_vec[1:num_QAOA_layers]
        betas = @view gammas_betas_combined_vec[num_QAOA_layers+1:end]

        proxy_amplitudes = QAOA_proxy(
            proxy,
            gammas,
            betas,
        )

        final_proxy_amplitudes = @view proxy_amplitudes[:,end]

        expectation = QAOA_proxy_expectation(
            proxy,
            final_proxy_amplitudes
        )
        current_time = time()

        if !isnothing(expectations)
            push!(expectations, (current_time, expectation))
        end

        return -expectation
    end

    return inverse_objective
end
