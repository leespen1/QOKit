"""
This is a modified version of QAOA_proxy.jl, which has additional 
arguments for h_peak (Float64) and distance (Int). 
This is for calling the QAOA proxy with different parameters for the proxy distribution we use.
    
Recommended reasonable settings for these parameters are: 
- h_peak = 1 << (num_qubits - 4) (this is the default)
- distance = num_qubits // 2 (this is the default)

An improvement for this would be to be able to pass in a distribution function instead of h_peak and distance, 
but this is not implemented here."""


"""
N(c'; d, c) from paper but instead of a multinomial distribution, we just approximate by a prism whose cross-sections at fixed distances are triangles
TODO: This only works for prob_edge = 0.5
"""
function number_of_costs_at_distance_proxy(cost_1::Int, cost_2::Int,
    distance::Int, num_constraints::Int, num_qubits::Int, h_peak::Float64)::Float64

    # Want distance to be between 0 and num_qubits//2 since further distance corresponds to being near the bitwise complement (which has the same cost)
    reflected_distance = distance
    if distance > div(num_qubits, 2)
    reflected_distance = num_qubits - distance
    end

    # Take the peak height at reflected_distance to be on the straight line between (0 or num_qubits, 1) and (num_qubits/2, h_peak)
    h_at_cost_2 = line_between(reflected_distance, 0, 1, num_qubits / 2, h_peak)
    # Let the peak height at reflected_distance occur where cost_2 is on the stright line between cost_1 and num_constraints/2
    center = line_between(reflected_distance, 0, cost_1, num_qubits / 2, num_constraints / 2)
    left = center - reflected_distance - 1
    right = center + reflected_distance + 1

    return triangle_value(cost_2, left, right, h_at_cost_2)
end

"""
Computes the sum inside the for loop of Algorithm 1 in paper using dumb approximations
"""
function compute_amplitude_sum(prev_amplitudes::AbstractVector{ComplexF64},
    gamma::Real, beta::Real, cost_1::Int, num_constraints::Int,
    num_qubits::Int, h_peak::Float64, distance::Int)::ComplexF64

    sum = 0
    for cost_2 in 0:num_constraints
    # Only use the provided distance
    beta_factor = (cos(beta) ^ (num_qubits - distance)) * ((-1im * sin(beta)) ^ distance)
    gamma_factor = exp(-1im * gamma * cost_2)
    num_costs_at_distance = number_of_costs_at_distance_proxy(cost_1, cost_2, distance, num_constraints, num_qubits, h_peak)
    sum += beta_factor * gamma_factor * prev_amplitudes[1+cost_2] * num_costs_at_distance
    end
    return sum
end

"""
Currently only implemented for prob_edge = 0.5
Now takes h_peak (Float64) and distance (Int) as arguments, which are used in number_of_costs_at_distance_proxy.
"""
function QAOA_proxy_peak_distance(p::Int, gamma::Vector{Float64}, beta::Vector{Float64}, num_constraints::Int, num_qubits::Int, h_peak::Float64, distance::Int, terms_to_drop_in_expectation::Int = 0)
    num_costs = num_constraints + 1
    amplitude_proxies = zeros(ComplexF64, p + 1, num_costs)
    init_amplitude = sqrt(1 / (1 << num_qubits))
    amplitude_proxies[1,:] .= init_amplitude # Memory inefficient, would be better to fill a column than a row
    
    for current_depth in 1:p
    for cost_1 in 0:num_costs-1
        amplitude_proxies[1+current_depth,1+cost_1] = compute_amplitude_sum(
        amplitude_proxies[current_depth,:], gamma[current_depth], beta[current_depth], cost_1, num_constraints, num_qubits, h_peak, distance
        )
    end
    end

    expected_proxy = 0
    for cost in terms_to_drop_in_expectation:num_costs-1
    expected_proxy += number_with_cost_proxy(cost, num_constraints, num_qubits) * (abs(amplitude_proxies[end,1+cost]) ^ 2) * cost
    end

    return amplitude_proxies, expected_proxy
end

"""
Convert numpy arrays to julia arrays before doing QAOA_proxy_peak_distance.
(should check whether this impacts performance)
"""
function QAOA_proxy_peak_distance(p, gamma::PyArray, beta::PyArray, num_constraints, num_qubits, h_peak, distance, terms_to_drop_in_expectation)
    return QAOA_proxy_peak_distance(
    pyconvert(Int, p),
    pyconvert(Vector, gamma),
    pyconvert(Vector, beta),
    pyconvert(Int, num_constraints),
    pyconvert(Int, num_qubits),
    pyconvert(Float64, h_peak),
    pyconvert(Int, distance),
    pyconvert(Int, terms_to_drop_in_expectation),
    )
end
