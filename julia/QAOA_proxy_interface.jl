"""
Computes the sum inside the for loop of Algorithm 1 in paper using dumb approximations
"""
function compute_amplitude_sum(
        proxy,
        prev_amplitudes::AbstractVector{ComplexF64},
        gamma::Real,
        beta::Real,
        cost_1::Integer
    )::ComplexF64

    @assert length(prev_amplitudes) == proxy.num_constraints + 1

    sum::ComplexF64 = 0
    # Changed loop order for efficiency. Check for same result
    for distance in 0:proxy.num_qubits 
        sinbeta, cosbeta = sincos(beta)
        beta_factor = (cosbeta ^ (proxy.num_qubits - distance)) * ((-1im * sinbeta) ^ distance)
        for cost_2 in 0:proxy.num_constraints
            gamma_factor = exp(-1im * gamma * cost_2)
            num_costs_at_distance = N_cost_distance_distribution(
                proxy, cost_1, distance, cost_2
            )
            sum += beta_factor * gamma_factor * prev_amplitudes[1+cost_2] * num_costs_at_distance
        end
    end
    return sum
end

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
(should check whether this impacts performance)
"""
function QAOA_proxy(proxy, gamma::PyArray, beta::PyArray)::Matrix{ComplexF64}
    return QAOA_proxy(
        proxy,
        pyconvert(Vector, gamma),
        pyconvert(Vector, beta),
    )
end

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
