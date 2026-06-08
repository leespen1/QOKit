################################################################################
#
# Julia implementation of real distribution functions from real_distribution.py
#
# These functions compute n(x; d, c) - the count of bitstrings at Hamming
# distance d with cost c from bitstring x - and N(c'; d, c) - the cost-averaged
# homogeneous distribution.
#
# This module only implements the functions that operate on costs arrays,
# NOT functions that depend on graph representations (NetworkX/QOKit).
# The costs arrays should be obtained from Python and passed to these functions.
#
# GPU ACCELERATION NOTES:
# The main computational bottleneck is get_real_distribution_from_costs which
# has O(2^n × 2^n) complexity. This can be GPU-accelerated by:
# 1. Using gpu_get_real_distribution_from_costs for CUDA arrays
# 2. The Hamming distance computation is trivially parallelizable
# 3. Atomic operations handle the histogram accumulation
#
################################################################################

export hamming_distance
export get_real_distribution_from_costs, get_homogeneous_distribution_from_costs
export get_homogeneous_distribution_from_costs_direct
export pad_to_shape, pad_to_match, pad_and_stack
export average_distributions, stddev_distributions, distributions_mean_and_stddev
export distribution_array_to_dict, get_pearson_correlation_coefficients


#==============================================================================#
#                          Core Distribution Functions                          #
#==============================================================================#

"""
    hamming_distance(x::Integer, y::Integer) -> Int

Compute the Hamming distance between two bitstrings represented as integers.

The Hamming distance is the number of positions at which the corresponding bits
differ. Uses Julia's built-in `count_ones` for efficiency.

# Examples
```julia
hamming_distance(0b000, 0b111)  # returns 3
hamming_distance(0b101, 0b111)  # returns 1
hamming_distance(0b010, 0b101)  # returns 3
```
"""
@inline function hamming_distance(x::Integer, y::Integer)::Int
    return count_ones(x ⊻ y)  # ⊻ is XOR, count_ones counts set bits
end


"""
    get_real_distribution_from_costs(costs, num_edges, num_vertices; max_num_edges=0) -> Array{Float64, 3}

Compute the real distribution n(x; d, c) from an array of costs.

n(x; d, c) gives the number of bitstrings y that have cost c and are at
Hamming distance d from bitstring x.

# Arguments
- `costs`: 1D array of costs for each bitstring (length 2^num_vertices)
- `num_edges`: Number of edges in the graph
- `num_vertices`: Number of vertices (qubits) in the graph
- `max_num_edges=0`: Optional, allocate extra space in cost axis for compatibility

# Returns
- 3D array of shape (num_bitstrings, num_distances, num_costs) where:
  - Index 1: bitstring x (0 to 2^num_vertices - 1)
  - Index 2: Hamming distance d (0 to num_vertices)
  - Index 3: cost c (0 to num_edges, or max_num_edges if larger)

# Notes
This is the main computational bottleneck with O(2^(2n)) complexity.
For GPU acceleration, use `gpu_get_real_distribution_from_costs`.
"""
function get_real_distribution_from_costs(
    costs::AbstractVector{<:Real},
    num_edges::Integer,
    num_vertices::Integer;
    max_num_edges::Integer=0
)::Array{Float64, 3}

    num_bitstrings = 2^num_vertices
    num_distances = num_vertices + 1
    num_costs = num_edges + 1
    cost_axis_size = max(num_costs, max_num_edges + 1)

    @assert length(costs) == num_bitstrings "Length of costs must equal 2^num_vertices"

    # Allocate output array (1-indexed in Julia, so size is exact)
    # Shape: (bitstring x, distance d, cost c)
    n_distribution = zeros(Float64, num_bitstrings, num_distances, cost_axis_size)

    # Main computation loop - O(2^(2n))
    # Each x writes to its own row n_distribution[x+1, :, :], so no contention
    @inbounds @threads for x in 0:(num_bitstrings - 1)
        for y in 0:(num_bitstrings - 1)
            d = hamming_distance(x, y)
            cost_y = Int(costs[y + 1])  # +1 for Julia 1-indexing

            # +1 for Julia 1-indexing on all array dimensions
            n_distribution[x + 1, d + 1, cost_y + 1] += 1.0
        end
    end

    return n_distribution
end


"""
    get_homogeneous_distribution_from_costs(costs, real_distribution; max_num_edges=0) -> Array{Float64, 3}

Compute the homogeneous distribution N(c'; d, c) by averaging n(x; d, c) over
all bitstrings x with cost c'.

# Arguments
- `costs`: 1D array of costs for each bitstring
- `real_distribution`: 3D array n(x; d, c) from `get_real_distribution_from_costs`
- `max_num_edges=0`: Optional, allocate extra space for compatibility

# Returns
- 3D array of shape (num_costs, num_distances, num_costs) representing N(c'; d, c)
"""
function get_homogeneous_distribution_from_costs(
    costs::AbstractVector{<:Real},
    real_distribution::AbstractArray{<:Real, 3};
    max_num_edges::Integer=0
)::Array{Float64, 3}

    num_bitstrings, num_distances, num_costs = size(real_distribution)
    cost_axis_size = max(num_costs, max_num_edges + 1)

    # Count cost occurrences for normalization (single-threaded, cheap)
    num_cost_occurrences = zeros(Int, num_costs)
    @inbounds for bitstring_i in 1:num_bitstrings
        bitstring_cost = Int(costs[bitstring_i]) + 1  # +1 for 1-indexing
        num_cost_occurrences[bitstring_cost] += 1
    end

    # Allocate one accumulator per thread to avoid write contention
    # (multiple bitstrings can share the same cost). Size by maxthreadid() so
    # the threadid()-indexed lookup is valid for any thread the loop runs on,
    # including interactive-pool threads whose ids exceed nthreads().
    nt = maxthreadid()
    thread_accumulators = [zeros(Float64, cost_axis_size, num_distances, cost_axis_size) for _ in 1:nt]

    # Sum n(x; d, c) over all x with same cost c', using thread-local accumulators
    @inbounds @threads :static for bitstring_i in 1:num_bitstrings
        tid = threadid()
        local_acc = thread_accumulators[tid]
        bitstring_cost = Int(costs[bitstring_i]) + 1  # +1 for 1-indexing

        for d in 1:num_distances
            for c in 1:num_costs
                local_acc[bitstring_cost, d, c] += real_distribution[bitstring_i, d, c]
            end
        end
    end

    # Reduce: sum all thread-local accumulators
    homogeneous_distribution = thread_accumulators[1]
    for t in 2:nt
        homogeneous_distribution .+= thread_accumulators[t]
    end

    # Normalize by number of bitstrings with each cost
    @inbounds for cost_idx in 1:num_costs
        if num_cost_occurrences[cost_idx] > 0
            homogeneous_distribution[cost_idx, :, :] ./= num_cost_occurrences[cost_idx]
        end
    end

    return homogeneous_distribution
end


"""
    get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices; max_num_edges=0) -> Array{Float64, 3}

Compute N(c'; d, c) directly from costs without storing the full n(x; d, c).

This is more memory-efficient for large systems as it avoids allocating the
O(2^n × n × m) real distribution array. Instead, it computes N(c'; d, c)
directly with O(m × n × m) storage.

# Arguments
- `costs`: 1D array of costs for each bitstring (length 2^num_vertices)
- `num_edges`: Number of edges in the graph
- `num_vertices`: Number of vertices in the graph
- `max_num_edges=0`: Optional, allocate extra space for compatibility

# Returns
- 3D array of shape (num_costs, num_distances, num_costs) representing N(c'; d, c)

# Notes
This function is the primary target for GPU acceleration since it combines
the two-step process into one and is more memory efficient.
"""
function get_homogeneous_distribution_from_costs_direct(
    costs::AbstractVector{<:Real},
    num_edges::Integer,
    num_vertices::Integer;
    max_num_edges::Integer=0
)::Array{Float64, 3}

    num_bitstrings = 2^num_vertices
    num_distances = num_vertices + 1
    num_costs = num_edges + 1
    cost_axis_size = max(num_costs, max_num_edges + 1)

    @assert length(costs) == num_bitstrings "Length of costs must equal 2^num_vertices"

    # Count occurrences of each cost (for normalization, single-threaded, cheap)
    num_cost_occurrences = zeros(Int, cost_axis_size)
    @inbounds for bitstring_i in 1:num_bitstrings
        cost_x = Int(costs[bitstring_i]) + 1
        num_cost_occurrences[cost_x] += 1
    end

    # Allocate one accumulator per thread to avoid write contention
    # (multiple x values can share the same cost_x). Size by maxthreadid() so
    # the threadid()-indexed lookup is valid for any thread the loop runs on,
    # including interactive-pool threads whose ids exceed nthreads().
    nt = maxthreadid()
    thread_accumulators = [zeros(Float64, cost_axis_size, num_distances, cost_axis_size) for _ in 1:nt]

    # Main computation: for each pair (x, y), accumulate to N(c(x); d(x,y), c(y))
    @inbounds @threads :static for x in 0:(num_bitstrings - 1)
        tid = threadid()
        local_acc = thread_accumulators[tid]
        cost_x = Int(costs[x + 1]) + 1  # c'

        for y in 0:(num_bitstrings - 1)
            d = hamming_distance(x, y) + 1  # d, +1 for 1-indexing
            cost_y = Int(costs[y + 1]) + 1  # c, +1 for 1-indexing

            local_acc[cost_x, d, cost_y] += 1.0
        end
    end

    # Reduce: sum all thread-local accumulators
    homogeneous_distribution = thread_accumulators[1]
    for t in 2:nt
        homogeneous_distribution .+= thread_accumulators[t]
    end

    # Normalize by number of bitstrings with each cost c'
    @inbounds for cost_idx in 1:cost_axis_size
        if num_cost_occurrences[cost_idx] > 0
            homogeneous_distribution[cost_idx, :, :] ./= num_cost_occurrences[cost_idx]
        end
    end

    return homogeneous_distribution
end


#==============================================================================#
#                          Utility Functions                                    #
#==============================================================================#

"""
    pad_to_shape(arr, target_shape) -> Array

Pad array with zeros to match target_shape.
"""
function pad_to_shape(arr::AbstractArray{T, N}, target_shape::NTuple{N, Int}) where {T, N}
    @assert all(size(arr) .<= target_shape) "arr shape must be <= target_shape elementwise"

    result = zeros(T, target_shape)
    indices = CartesianIndices(size(arr))
    result[indices] .= arr

    return result
end


"""
    pad_to_match(a, b) -> (Array, Array)

Pad the smaller array to match the shape of the larger array.
Returns both arrays with matching shapes.
"""
function pad_to_match(a::AbstractArray{Ta, N}, b::AbstractArray{Tb, N}) where {Ta, Tb, N}
    shape_a = size(a)
    shape_b = size(b)

    if shape_a == shape_b
        return a, b
    end

    if all(shape_a .<= shape_b)
        return pad_to_shape(a, shape_b), b
    elseif all(shape_b .<= shape_a)
        return a, pad_to_shape(b, shape_a)
    else
        throw(ArgumentError("Shapes are not broadcast-compatible for padding"))
    end
end


"""
    pad_and_stack(arrays) -> Array{T, 4}

Pad all 3D arrays to the largest shape among them, then stack along a new first axis.
Returns array of shape (num_arrays, dim1, dim2, dim3).
"""
function pad_and_stack(arrays::AbstractVector{<:AbstractArray{T, 3}}) where T
    @assert all(ndims(arr) == 3 for arr in arrays) "All arrays must be 3D"

    max_shape = (
        maximum(size(arr, 1) for arr in arrays),
        maximum(size(arr, 2) for arr in arrays),
        maximum(size(arr, 3) for arr in arrays)
    )

    padded = [pad_to_shape(arr, max_shape) for arr in arrays]
    # Stack along first dimension: result is (num_arrays, dim1, dim2, dim3)
    return stack(padded; dims=1)
end


"""
    average_distributions(distributions) -> Array{Float64, 3}

Compute the element-wise average of multiple 3D distributions.
Arrays are padded to match the largest shape.
"""
function average_distributions(distributions::AbstractVector{<:AbstractArray{<:Real, 3}})
    stacked = pad_and_stack(distributions)
    return dropdims(mean(stacked, dims=1), dims=1)
end


"""
    stddev_distributions(distributions) -> Array{Float64, 3}

Compute the element-wise standard deviation of multiple 3D distributions.
Uses population std (N denominator) to match NumPy's default behavior.
"""
function stddev_distributions(distributions::AbstractVector{<:AbstractArray{<:Real, 3}})
    stacked = pad_and_stack(distributions)
    return dropdims(std(stacked; dims=1, corrected=false), dims=1)
end


"""
    distributions_mean_and_stddev(distributions) -> (mean, stddev)

Compute both mean and standard deviation of multiple 3D distributions.
Returns tuple of (mean_array, stddev_array).
Uses population std (N denominator) to match NumPy's default behavior.
"""
function distributions_mean_and_stddev(distributions::AbstractVector{<:AbstractArray{<:Real, 3}})
    stacked = pad_and_stack(distributions)
    m = dropdims(mean(stacked, dims=1), dims=1)
    s = dropdims(std(stacked; dims=1, corrected=false), dims=1)
    return m, s
end


"""
    distribution_array_to_dict(distribution_array) -> Dict

Convert a 3D distribution array to a dictionary where keys are (i, j, k) tuples
and values are the non-zero distribution values.

Useful for sparse representations and cross-language interoperability.
"""
function distribution_array_to_dict(distribution_array::AbstractArray{T, 3}) where T
    result = Dict{Tuple{Int, Int, Int}, T}()

    for idx in CartesianIndices(distribution_array)
        val = distribution_array[idx]
        if val != zero(T)
            # Convert to 0-indexed tuple for Python compatibility
            result[(idx[1] - 1, idx[2] - 1, idx[3] - 1)] = val
        end
    end

    return result
end


"""
    get_pearson_correlation_coefficients(homodist1, homodist2) -> Vector{Float64}

Compute Pearson correlation coefficients between corresponding 2D slices
of two homogeneous distributions.

For each cost c', computes the correlation between homodist1[c', :, :] and
homodist2[c', :, :] (flattened to 1D vectors).
"""
function get_pearson_correlation_coefficients(
    homodist1::AbstractArray{<:Real, 3},
    homodist2::AbstractArray{<:Real, 3}
)
    homodist1, homodist2 = pad_to_match(homodist1, homodist2)

    function pearson_corr(arr1, arr2)
        v1 = vec(arr1)
        v2 = vec(arr2)
        # cor computes Pearson correlation
        return cor(v1, v2)
    end

    num_costs = size(homodist1, 1)
    return [pearson_corr(homodist1[i, :, :], homodist2[i, :, :]) for i in 1:num_costs]
end
