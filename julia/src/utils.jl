abstract type AbstractProxy end

"""
Allocate an empty array full of undefined elements, to store the homogeneous
distribution in.

Note that order is (c', c, d), NOT (c', d, c).
"""
function allocate_homodist(proxy::AbstractProxy)
    distribution = allocate_homodist(proxy.num_constraints, proxy.num_qubits)
    return distribution
end

function allocate_homodist(num_constraints::Integer, num_qubits::Integer)
    num_costs = num_constraints + 1
    num_distances = num_qubits + 1
    distribution = Array{Float64, 3}(undef, num_costs, num_distances, num_costs)
    return distribution
end

"""
Construct N(c',d,c) from proxy.
"""
function cpu_compute_homodist(proxy)
    # These three are 'vectors', but along different dimensions, for broadcasting to work
    costs_prime = collect(0:proxy.num_constraints)
    distances = reshape(collect(0:proxy.num_qubits), 1, :)
    costs_unprime = reshape(costs_prime, 1, 1, :)
    N(c_prime, c, d) = N_cost_distance_distribution(proxy, c_prime, d, c)
    homodist = N.(costs_prime, costs_unprime, distances)
    return homodist
end


function gpu_compute_homodist(proxy)
    # These three are 'vectors', but along different dimensions, for broadcasting to work
    costs_prime = collect(0:proxy.num_constraints) |> CuArray
    distances = reshape(collect(0:proxy.num_qubits), 1, :) |> CuArray
    costs_unprime = reshape(costs_prime, 1, 1, :) |> CuArray
    N(c_prime, c, d) = N_cost_distance_distribution(proxy, c_prime, d, c)
    homodist = N.(costs_prime, costs_unprime, distances)
    return homodist
end

"""
- `proxies`: A vector of proxies 
- `batch_size`: number of proxies to do in one kernel call to GPU (limit as
needed to meet memory reguirements).

Assumed that all proxies have some num_qubits and num_constraints

Consier lowering precision in the future. I think 32-bit should be fine.
"""
function gpu_multi_proxy_mse(
    proxies::AbstractVector{<: AbstractProxy},
    sampled_homodist::AbstractArray{<: Real, 3};
    batch_size::Integer=length(proxies),
    normalize=true,
    show_progress=false
)
    proxy = proxies[1] 
     
    # Could make costs arrays of UInt16's to save space.
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
        if normalize # Normalize each homodist
            homodist_gpu_volumes = sum(sampled_homodist_gpu, dims=(1,2,3))
            homodists_gpu ./= homodist_gpu_volumes # Vectorized, should divide each 3D slice
        end
        homodists_gpu .-= sampled_homodist_gpu # sampled_homodist is 3D, this will be repeated over 4th dimension

        mse_batch_gpu = mapreduce(x -> x*x, +, homodists_gpu, dims=(1,2,3)) # Reduce first 3 dims, leave 4th dim
        mse_batch_gpu ./= num_elements_in_homodist

        push!(mse_batches, Array(vec(mse_batch_gpu)))
    end
    return mse_vec = mse_batches |> Iterators.flatten |> collect
end

function cpu_multi_proxy_mse(
    proxies::AbstractVector{<: AbstractProxy},
    sampled_homodist::AbstractArray{<: Real, 3};
    batch_size::Integer=length(proxies),
    normalize=true,
    show_progress=false,
)
    if normalize
        sampled_homodist_volume = sum(sampled_homodist_gpu)
        sampled_homodist_gpu ./= sampled_homodist_volume
    end
     
    # Could make costs arrays of UInt16's.
    proxy = proxies[1] 
    costs_prime = 0:proxy.num_constraints
    distances = reshape(0:proxy.num_qubits, (1, :))
    costs_unprime = reshape(costs_prime, (1, 1, :))
    num_elements_in_homodist = (1+proxy.num_qubits)*(1+proxy.num_constraints)^2
    mse_batches = Vector{Float64}[]
    @showprogress enabled=show_progress for (i, proxy_batch) in enumerate(Iterators.partition(proxies, batch_size))
        proxy_batch_reshaped = reshape(proxy_batch, (1,1,1,:)) 
        homodists = N_cost_distance_distribution.(proxy_batch_reshaped, costs_prime, distances, costs_unprime)  
        if normalize # Normalize each homodist
            homodist_volumes = sum(homodists, dims=(1,2,3))
            homodists ./= homodist_volumes # Vectorized, should divide each 3D slice
        end
        homodists .-= sampled_homodist # sampled_homodist is 3D, this will be repeated over 4th dimension
        mse_batch = mapreduce(x -> x*x, +, homodists, dims=(1,2,3)) # Reduce first 3 dims, leave 4th dim
        mse_batch ./= num_elements_in_homodist

        push!(mse_batches, vec(mse_batch))
    end
    return mse_vec = mse_batches |> Iterators.flatten |> collect
end

"""
Compute the sum of the squared errors between two arrays, given by
∑ᵢ (aᵢ- bᵢ)² 
To geat the mean squared error, simply divide by the number of elements.
"""
@inline function sum_squared_error(array1, array2)
    return mapreduce((x,y) -> (x-y)^2, +, array1, array2)
end


