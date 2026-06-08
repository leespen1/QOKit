module JuliaQAOACUDAExt

using CUDA
using JuliaQAOA
using JuliaQAOA: AbstractProxy, N_cost_distance_distribution, _expand
using ProgressMeter: @showprogress


# ─── gpu_compute_homodist ─────────────────────────────────────────────────────

function JuliaQAOA.gpu_compute_homodist(proxy)
    costs_prime = collect(0:proxy.num_constraints) |> CuArray
    distances = reshape(collect(0:proxy.num_qubits), 1, :) |> CuArray
    costs_unprime = reshape(costs_prime, 1, 1, :) |> CuArray
    N(c_prime, d, c) = N_cost_distance_distribution(proxy, c_prime, d, c)
    homodist = N.(costs_prime, distances, costs_unprime)
    return homodist
end


# ─── gpu_multi_proxy_mse ──────────────────────────────────────────────────────

function JuliaQAOA.gpu_multi_proxy_mse(
    proxies::AbstractVector{<: AbstractProxy},
    sampled_homodist::AbstractArray{<: Real, 3};
    batch_size::Integer=length(proxies),
    normalize=true,
    show_progress=false
)
    proxy = proxies[1]

    costs_prime = collect(0:proxy.num_constraints) |> CuArray
    distances = reshape(collect(0:proxy.num_qubits), (1, :)) |> CuArray
    costs_unprime = reshape(costs_prime, (1, 1, :)) |> CuArray
    num_elements_in_homodist = (1+proxy.num_qubits)*(1+proxy.num_constraints)^2
    mse_batches = Vector{Float64}[]
    sampled_homodist_gpu = CuArray(sampled_homodist)

    if normalize
        sampled_homodist_volume = sum(sampled_homodist_gpu)
        sampled_homodist_gpu ./= sampled_homodist_volume
    end

    @showprogress enabled=show_progress for (i, proxy_batch_gpu) in enumerate(Iterators.partition(proxies, batch_size))
        proxy_batch_gpu_reshaped = reshape(proxy_batch_gpu, (1,1,1,:)) |> CuArray

        homodists_gpu = N_cost_distance_distribution.(
            proxy_batch_gpu_reshaped, costs_prime, distances, costs_unprime
        )

        if normalize
            homodist_gpu_volumes = sum(homodists_gpu, dims=(1,2,3))
            homodists_gpu ./= homodist_gpu_volumes
        end

        homodists_gpu .-= sampled_homodist_gpu
        mse_batch_gpu = mapreduce(x -> x*x, +, homodists_gpu, dims=(1,2,3))
        mse_batch_gpu ./= num_elements_in_homodist

        push!(mse_batches, Array(vec(mse_batch_gpu)))
    end
    return mse_vec = mse_batches |> Iterators.flatten |> collect
end


# ─── GPU cost distributions (CUDA-specific) ──────────────────────────────────

"""
    gpu_get_real_distribution_from_costs(costs, num_edges, num_vertices; max_num_edges=0) -> CuArray{Float64, 3}

GPU-accelerated version of get_real_distribution_from_costs.

Uses CUDA kernels to parallelize the O(2^(2n)) computation across GPU threads.
Each thread handles one (x, y) pair, using atomic operations to accumulate counts.

# Arguments
Same as `get_real_distribution_from_costs`

# Returns
- CuArray of shape (num_bitstrings, num_distances, num_costs)
"""
function JuliaQAOA.gpu_get_real_distribution_from_costs(
    costs::AbstractVector{<:Real},
    num_edges::Integer,
    num_vertices::Integer;
    max_num_edges::Integer=0
)::CuArray{Float64, 3}

    num_bitstrings = 2^num_vertices
    num_distances = num_vertices + 1
    num_costs = num_edges + 1
    cost_axis_size = max(num_costs, max_num_edges + 1)

    @assert length(costs) == num_bitstrings "Length of costs must equal 2^num_vertices"

    costs_gpu = CuArray{Int32}(Int32.(costs))
    n_distribution = CUDA.zeros(Float64, num_bitstrings, num_distances, cost_axis_size)

    threads_per_block = 256
    num_pairs = num_bitstrings * num_bitstrings
    blocks = cld(num_pairs, threads_per_block)

    @cuda threads=threads_per_block blocks=blocks _real_distribution_kernel!(
        n_distribution, costs_gpu, num_bitstrings, num_vertices
    )

    return n_distribution
end

"""
CUDA kernel for computing real distribution n(x; d, c).
Each thread handles one (x, y) pair.
"""
function _real_distribution_kernel!(
    n_distribution::CuDeviceArray{Float64, 3},
    costs::CuDeviceArray{Int32, 1},
    num_bitstrings::Int,
    num_vertices::Int
)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    num_pairs = num_bitstrings * num_bitstrings

    if idx <= num_pairs
        x = (idx - 1) ÷ num_bitstrings
        y = (idx - 1) % num_bitstrings
        d = count_ones(x ⊻ y)
        cost_y = costs[y + 1]
        CUDA.@atomic n_distribution[x + 1, d + 1, cost_y + 1] += 1.0
    end

    return nothing
end


"""
    gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices; max_num_edges=0) -> CuArray{Float64, 3}

GPU-accelerated version that computes N(c'; d, c) directly.

This is more memory-efficient than computing the full n(x; d, c) first,
as it only allocates O(m × n × m) storage instead of O(2^n × n × m).

Uses two optimizations over a naive GPU implementation:
1. Int32 atomics instead of Float64 (faster native support and higher throughput)
2. Privatized global memory: multiple copies of the output array reduce
   atomic contention by spreading thread blocks across independent copies.

# Arguments
Same as `get_homogeneous_distribution_from_costs_direct`

# Returns
- CuArray of shape (num_costs, num_distances, num_costs)
"""
function JuliaQAOA.gpu_get_homogeneous_distribution_from_costs_direct(
    costs::AbstractVector{<:Real},
    num_edges::Integer,
    num_vertices::Integer;
    max_num_edges::Integer=0
)::CuArray{Float64, 3}

    num_bitstrings = 2^num_vertices
    num_distances = num_vertices + 1
    num_costs = num_edges + 1
    cost_axis_size = max(num_costs, max_num_edges + 1)

    @assert length(costs) == num_bitstrings "Length of costs must equal 2^num_vertices"

    costs_gpu = CuArray{Int32}(Int32.(costs))

    dev = CUDA.device()
    num_sms = CUDA.attribute(dev, CUDA.DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT)
    num_copies = max(min(num_sms, 64), 1)

    privatized_counts = CUDA.zeros(Int32, cost_axis_size, num_distances, cost_axis_size, num_copies)

    num_cost_occurrences = zeros(Float64, cost_axis_size)
    @inbounds for bitstring_i in 1:num_bitstrings
        cost_x = Int(costs[bitstring_i]) + 1
        num_cost_occurrences[cost_x] += 1.0
    end

    threads_per_block = 256
    num_pairs = num_bitstrings * num_bitstrings
    blocks = cld(num_pairs, threads_per_block)

    @cuda threads=threads_per_block blocks=blocks _homogeneous_distribution_kernel_priv!(
        privatized_counts, costs_gpu, num_bitstrings, num_vertices, num_copies
    )

    CUDA.synchronize()

    total_counts = dropdims(sum(privatized_counts, dims=4), dims=4)
    homogeneous_distribution = Float64.(total_counts)

    num_cost_occurrences_gpu = CuArray(num_cost_occurrences)
    normalizers = reshape(num_cost_occurrences_gpu, :, 1, 1)
    normalizers = max.(normalizers, 1.0)
    homogeneous_distribution ./= normalizers

    return homogeneous_distribution
end

"""
CUDA kernel for computing homogeneous distribution N(c'; d, c) directly,
using privatized global memory to reduce atomic contention.
Each thread handles one (x, y) pair and writes to its block's assigned copy.
Uses Int32 atomics for maximum throughput.
"""
function _homogeneous_distribution_kernel_priv!(
    privatized_counts::CuDeviceArray{Int32, 4},
    costs::CuDeviceArray{Int32, 1},
    num_bitstrings::Int,
    num_vertices::Int,
    num_copies::Int
)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    num_pairs = num_bitstrings * num_bitstrings

    if idx <= num_pairs
        x = (idx - 1) ÷ num_bitstrings
        y = (idx - 1) % num_bitstrings
        cost_x = costs[x + 1]
        cost_y = costs[y + 1]
        d = count_ones(x ⊻ y)
        copy_idx = ((blockIdx().x - 1) % num_copies) + 1
        CUDA.@atomic privatized_counts[cost_x + 1, d + 1, cost_y + 1, copy_idx] += Int32(1)
    end

    return nothing
end


# ─── Warp-shuffle X-mixer ────────────────────────────────────────────────────

include("batched_furx_cuda.jl")


end # module
