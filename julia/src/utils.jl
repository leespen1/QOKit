abstract type AbstractProxy end

"""
Allocate an empty array full of undefined elements, to store the homogeneous
distribution in.
"""
function allocate_homodist(proxy::AbstractProxy)
    num_constraints = proxy.num_constraints
    num_costs = num_constraints + 1
    num_distances = proxy.num_qubits + 1
    distribution = Array{Float64, 3}(undef, num_costs, num_costs, num_distances)
    return distribution
end

"""
Construct N(c',c,d) from proxy.

NOTE THE DIFFERENT ORDERING THAN USUAL! IT IS NOT N(c',d,c)!
"""
function efficient_order_homodist(
    proxy,
    n_threads = Base.Threads.nthreads()
)
    num_constraints = proxy.num_constraints
    num_costs = num_constraints + 1
    num_distances = proxy.num_qubits + 1
    distribution = zeros(Float64, num_costs, num_costs, num_distances)
    # Could make this parallel (on GPU too, with easy with tensor contractions)
    # But first make it work and profile. Don't optimize prematurely
    n_elements = length(distribution)
    #n_threads = 8
    chunk_size = div(n_elements, n_threads, RoundUp)
    cartesian_indices = CartesianIndices(distribution)
    @threads for t in 1:n_threads
        start_id = 1 + (t-1)*(chunk_size)
        end_id = min(t*chunk_size, n_elements)
        for id in start_id:end_id
            cart_id = cartesian_indices[id] # Use Cartesian indices to get integer costs and distance
            cost_prime, cost, distance = Tuple(cart_id) .- 1
            distribution[id] = N_cost_distance_distribution(proxy, cost_prime, distance, cost)
        end
    end
    return distribution
end

"""
Obtained very similar performance compared to cartesian indices approach,
only 2% faster
"""
function simple_homodist(proxy)
    # These three are 'vectors', but along different dimensions, for broadcasting to work
    costs_prime = collect(0:proxy.num_constraints)
    costs_unprime = reshape(costs_prime, 1, :)
    distances = reshape(collect(0:proxy.num_qubits), 1, 1, :)
    N(c_prime, c, d) = N_cost_distance_distribution(proxy, c_prime, d, c)
    homodist = N.(costs_prime, costs_unprime, distances)
    return homodist
end


function gpu_simple_homodist(proxy)
    # These three are 'vectors', but along different dimensions, for broadcasting to work
    costs_prime = collect(0:proxy.num_constraints) |> CuArray
    costs_unprime = reshape(costs_prime, 1, :) |> CuArray
    distances = reshape(collect(0:proxy.num_qubits), 1, 1, :) |> CuArray
    N(c_prime, c, d) = N_cost_distance_distribution(proxy, c_prime, d, c)
    homodist = N.(costs_prime, costs_unprime, distances)
    return homodist
end

"""
- `proxies`: A vector of proxies 
- `batch_size`: number of proxies to do in one kernel call to GPU (limit as
needed to meet memory reguirements).

Assumed that all proxies have some num_qubits and num_constraints
"""
function gpu_triangle_proxy_sweep_homodist(
    proxies::AbstractVector{<: AbstractProxy},
    sampled_homodist::AbstractArray{<: Real, 3};
    batch_size::Integer=length(proxies)
)
     
    # Could make costs arrays of UInt16's.
    costs_prime = collect(0:proxy.num_constraints) |> CuArray
    costs_unprime = reshape(costs_prime, (1, :)) |> CuArray
    distances = reshape(collect(0:proxy.num_qubits), (1, 1, :)) |> CuArray
    num_elements_in_homodist = prod(length, (costs_prime, costs_unprime, distances))
    mse_batches = Vector{Float64}[]
    sampled_homodist_gpu = CuArray(sampled_homodist)
    for (i, proxy_batch) in enumerate(Iterators.partition(proxies, batch_size))
        proxy_batch_gpu = reshape(proxy_batch, (1,1,1,:)) |> CuArray
        homodists_gpu = N_cost_distance_distribution.(
            proxy_batch_gpu, costs_prime, distances, costs_unprime
        )  
        homodists_gpu .-= sampled_homodist_gpu # sampled_homodist is 3D, this will be repeated over 4th dimension
        mse_batch_gpu = mapreduce(x -> x*x, +, homodists, dims=(1,2,3)) # Reduce first 3 dims, leave 4th dim
        mse_batch_gpu ./= num_elements_in_homodist

        push!(mse_batches, Array(vec(mse_batch_gpu)))
    end
    return mse_vec = mse_batches |> Iterators.flatten |> collect
end

function cpu_triangle_proxy_sweep_homodist(
    proxies::AbstractVector{<: AbstractProxy},
    sampled_homodist::AbstractArray{<: Real, 3};
    batch_size::Integer=length(proxies)
)
     
    # Could make costs arrays of UInt16's.
    proxy = proxies[1] 
    costs_prime = 0:proxy.num_constraints
    costs_unprime = reshape(costs_prime, (1, :))
    distances = reshape(0:proxy.num_qubits, (1, 1, :))
    num_elements_in_homodist = prod(length, (costs_prime, costs_unprime, distances))
    mse_batches = Vector{Float64}[]
    for (i, proxy_batch) in enumerate(Iterators.partition(proxies, batch_size))
        proxy_batch_gpu = reshape(proxy_batch, (1,1,1,:)) 
        homodists = N_cost_distance_distribution.(proxy_batch_gpu, costs_prime, distances, costs_unprime)  
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


