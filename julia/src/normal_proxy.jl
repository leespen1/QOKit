struct NormalProxy
    num_constraints::Int64
    num_qubits::Int64
    cost_mean::Float64
    cov_1::Float64
    cov_2::Float64
end

function P_cost_distribution(proxy::NormalProxy, cost::Integer)::Float64
    prob_cost_mean = proxy.num_qubits / 2
    prob_cost_cov = proxy.num_qubits / 4
    normal_distribution = Normal(prob_cost_mean, prob_cost_cov)
    return pdf(normal_distribution, cost)
end

function N_cost_distribution(proxy::NormalProxy, cost::Integer)::Float64
    scale = 1 << proxy.num_qubits
    return P_cost_distribution(proxy,cost) * scale
end

# N(c'; d, c) from paper
function N_cost_distance_distribution(proxy::NormalProxy, cost_1::Integer, distance::Integer, cost_2::Integer)::Float64
    distance_mean = proxy.num_qubits / 2

    # Note: A = [a b;-b a], B = [c 0;0 d] implies A*B*A^(-1) is symmetric

    # Use static matrices to avoid dynamic memory allocations
    P = @SMatrix [(cost_1-proxy.cost_mean)  distance_mean;
                  -distance_mean            (cost_1 - proxy.cost_mean)]

    cov_diag_mat = @SMatrix [proxy.cov_1 0; 0 proxy.cov_2]
    cov_mat = (P * cov_diag_mat) / P # Right-hand solve, don't compute inverse directly
    cov_mat_sym = @SMatrix [cov_mat[1,1] cov_mat[1,2]; cov_mat[1,2] cov_mat[2,2]] # Enforce cov_mat is symmetric, avoid floating point erros
    means = SVector(proxy.cost_mean, distance_mean)
    multivariate_distribution = MvNormal(means, cov_mat_sym)
    scale = (1 << proxy.num_qubits)
    return scale * pdf(multivariate_distribution, SVector(cost_2, distance))
end
