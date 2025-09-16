"""
Construct N(c',c,d) from proxy.

NOTE THE DIFFERENT ORDERING THAN USUAL! IT IS NOT N(c',d,c)!
"""
function efficient_order_homodist(proxy)
    num_constraints = proxy.num_constraints
    num_costs = num_constraints + 1
    num_distances = proxy.num_qubits + 1
    distribution = zeros(Float64, num_costs, num_costs, num_distances)
    # Could make this parallel (on GPU too, with easy with tensor contractions)
    # But first make it work and profile. Don't optimize prematurely
    n_elements = length(distribution)
    n_threads = Base.Threads.nthreads()
    chunk_size = div(n_elements, n_threads, RoundUp)
    cartesian_indices = CartesianIndices(distribution)
    @threads for t in 1:n_threads
        println("Thread ", t, " of ", n_threads)
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
