abstract type AbstractProxy end

Base.iterate(proxy::AbstractProxy) = (proxy, nothing)
Base.iterate(proxy::AbstractProxy, state) = nothing
Base.length(proxy::AbstractProxy) = 1

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
    distribution = fill(1.0, num_costs, num_distances, num_costs)
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
    N(c_prime, d, c) = N_cost_distance_distribution(proxy, c_prime, d, c)
    homodist = N.(costs_prime, distances, costs_unprime)
    return homodist
end


function cpu_multi_proxy_mse(
    proxies::AbstractVector{<: AbstractProxy},
    sampled_homodist::AbstractArray{<: Real, 3};
    batch_size::Integer=length(proxies),
    normalize=true,
    show_progress=false,
)
    proxy = proxies[1] 
     
    # Could make costs arrays of UInt16's to save space.
    costs_prime = collect(0:proxy.num_constraints)
    distances = reshape(collect(0:proxy.num_qubits), (1, :))
    costs_unprime = reshape(costs_prime, (1, 1, :))
    num_elements_in_homodist = (1+proxy.num_qubits)*(1+proxy.num_constraints)^2
    mse_batches = Vector{Float64}[]

    if normalize
        sampled_homodist_volume = sum(sampled_homodist)
        sampled_homodist ./= sampled_homodist_volume
    end

    @showprogress enabled=show_progress for (i, proxy_batch) in enumerate(Iterators.partition(proxies, batch_size))
        proxy_batch_reshaped = reshape(proxy_batch, (1,1,1,:)) 

        homodists = N_cost_distance_distribution.(
            proxy_batch_reshaped, costs_prime, distances, costs_unprime
        )  

        if normalize # Normalize each homodist
            homodist_volumes = sum(homodists, dims=(1,2,3))
            homodists ./= homodist_volumes # Vectorized, should divide each 3D slice
        end

        # TODO mse_batch computation could be reduced to one mapreduce call
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


