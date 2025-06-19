struct HardCodedTriangleProxy
    num_constraints::Int64
    num_qubits::Int64
end

"""
Gives the y-value at x=current_time on the line between (start_time, start_value) and (end_time, end_value)
"""
function line_between(current_time::Real, start_time::Real,
        start_value::Real, end_time::Real, end_value::Real)::Float64

    # Goes from 0 to 1 as current_time goes from start_time to end_time
    relative_time = (current_time - start_time) / (end_time - start_time)

    # Goes from start_value to end_value as relative_time goes from 0 to 1
    return (1 - relative_time) * start_value + relative_time * end_value
end

raw"""
             /\height
            /  \
 _ _ _ left/    \right _ _ _ given x, returns the corresponding y value on the preceeding curve
"""
function triangle_value(x::Integer, left::Real, right::Real, height::Real)::Float64
    return max(0, min(x - left, right - x) * 2 * height / (right - left))
end

"""
P(c') from paper but dumber
"""
function P_cost_distribution(proxy::HardCodedTriangleProxy, cost::Integer)::Float64
    return 4 / ((proxy.num_constraints + 1) ^ 2) * min(cost + 1, proxy.num_constraints + 1 - cost)
end

"""
N(c') from paper but dumber
"""
function N_cost_distribution(proxy::HardCodedTriangleProxy, cost::Integer)::Float64
    scale = 1 << proxy.num_qubits
    return P_cost_distribution(proxy, cost) * scale
end

"""
N(c'; d, c) from paper but instead of a multinomial distribution, we just approximate by a prism whose cross-sections at fixed distances are triangles
TODO: This only works for prob_edge = 0.5
"""
function N_cost_distance_distribution(proxy::HardCodedTriangleProxy,
        cost_1::Integer, distance::Integer, cost_2::Integer)::Float64

    # Want distance to be between 0 and num_qubits//2 since further distance corresponds to being near the bitwise complement (which has the same cost)
    reflected_distance = distance
    if distance > div(proxy.num_qubits, 2)
        reflected_distance = proxy.num_qubits - distance
    end

    # Approximate the peak value of the paper's multinomial distribution (roughly)
    h_peak = 1 << (proxy.num_qubits - 4)
    # Take the peak height at reflected_distance to be on the straight line between (0 or num_qubits, 1) and (num_qubits/2, h_peak)
    h_at_cost_2 = line_between(reflected_distance, 0, 1, proxy.num_qubits / 2, h_peak)
    # Let the peak height at reflected_distance occur where cost_2 is on the stright line between cost_1 and num_constraints/2
    center = line_between(reflected_distance, 0, cost_1, proxy.num_qubits / 2, proxy.num_constraints / 2)
    left = center - reflected_distance - 1
    right = center + reflected_distance + 1

    return triangle_value(cost_2, left, right, h_at_cost_2)
end
