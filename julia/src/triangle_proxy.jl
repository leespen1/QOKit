#=
This file implements the parameterized "triangle" distribution, as well as the
original "hard-coded triangle" distribution.
=#

"""
Class for implementing the "Triangle" parameterized proxy for QAOA.

Required arguments:
- num_constraints: int
- num_qubits: int
Argument with defaults (optional):
- h_tweak_sub
- hc_tweak_add
- l_tweak_mul
- r_tweak_mul

For a graph with M edges and N vertices, we create a TriangleProxy like this:
    proxy = TriangleProxy(M, N, h_tweak_sub, hc_tweak_add, l_tweak_mul, r_tweak_mul)

"""
struct TriangleProxy
    num_constraints::Int64
    num_qubits::Int64
    h_tweak_sub::Float64  # Shifts the peak of the pyramid down (Default 0)
    hc_tweak_add::Float64 # Moves the cost_2 of the peak to the right (Default 0)
    l_tweak_mul::Float64  # Defines the (inverse of the) slope of the left side of the pyramid (Default 1)
    r_tweak_mul::Float64  # Defines the (inverse of the) slope of the right side of the pyramid (Default 1)
    function TriangleProxy( # For default arguements
            num_constraints,
            num_qubits,
            h_tweak_sub=0,
            hc_tweak_add=0,
            l_tweak_mul=1,
            r_tweak_mul=1
        )
        new(
            num_constraints,
            num_qubits,
            h_tweak_sub,
            hc_tweak_add,
            l_tweak_mul,
            r_tweak_mul
        )
    end
end

"""
P(c') from paper
"""
function P_cost_distribution(proxy::TriangleProxy, cost::Integer)::Float64
    return 4 / ((proxy.num_constraints + 1) ^ 2) * min(cost + 1, proxy.num_constraints + 1 - cost)
end

"""
N(c') from paper
"""
function N_cost_distribution(proxy::TriangleProxy, cost::Integer)::Float64
    scale = 1 << proxy.num_qubits
    return P_cost_distribution(proxy, cost) * scale
end


"""
N(c'; d, c) from paper
"""
function N_cost_distance_distribution(proxy::TriangleProxy, cost_1::Integer, distance::Integer, cost_2::Integer)::Float64
    # Want distance to be between 0 and proxy.num_qubits//2 since further distance corresponds to being near the bitwise complement (which has the same cost)
    reflected_distance = (distance > div(proxy.num_qubits, 2)) ? proxy.num_qubits - distance : distance

    # Approximate the peak value of the paper's multinomial distribution (roughly)
    h_peak = (1 << (proxy.num_qubits - 4)) - proxy.h_tweak_sub
    center_at_h_peak = proxy.num_constraints / 2 + proxy.hc_tweak_add
    # Take the peak height at reflected_distance to be on the straight line between (0 or proxy.num_qubits, 1) and (proxy.num_qubits/2, h_peak)
    h_at_cost_2 = line_between(reflected_distance, 0, 1, proxy.num_qubits / 2, h_peak)
    # Let the peak height at reflected_distance occur where cost_2 is on the stright line between cost_1 and proxy.num_constraints/2
    center = line_between(reflected_distance, 0, cost_1, proxy.num_qubits / 2, center_at_h_peak)
    left = center - proxy.l_tweak_mul * reflected_distance - 1
    right = center + proxy.r_tweak_mul * reflected_distance + 1

    return triangle_value(cost_2, left, right, h_at_cost_2)
end


"""
Type for the original, hard-coded triangle distribution.

This type is now obselete, because the same result can be achieved using
the parameterized triangle proxy with the default parameters.
"""
struct HardCodedTriangleProxy
    num_constraints::Int64
    num_qubits::Int64
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
N(c'; d, c) from paper but instead of a multinomial distribution, we just
approximate by a prism whose cross-sections at fixed distances are triangles
"""
function N_cost_distance_distribution(proxy::HardCodedTriangleProxy,
        cost_1::Integer, distance::Integer, cost_2::Integer)::Float64
    # Want distance to be between 0 and proxy.num_qubits//2 since further distance corresponds to being near the bitwise complement (which has the same cost)
    reflected_distance = (distance > div(proxy.num_qubits, 2)) ? proxy.num_qubits - distance : distance

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


"""
Gives the y-value at x=current_time on the line between (start_time,
start_value) and (end_time, end_value)
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
