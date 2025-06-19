using Distributions, StaticArrays, PythonCall

#"""
#This file implements the QAOA proxy algorithm for MaxCut from:
#https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.6.023171
#"""

struct PaperProxy
    num_constraints::Int64
    num_qubits::Int64
    prob_edge::Float64
    function PaperProxy(num_constraints, num_qubits, prob_edge=0.5)
        new(num_constraints, num_qubits, prob_edge)
    end
end



"""
P(c') from paper
"""
function P_cost_distribution(proxy::PaperProxy, cost::Integer)::Float64
    # Binomial distribution could be moved to time of proxy-creation,
    # but benchmarks suggest savings would be < 2%
    binomial_distribution = Binomial(proxy.num_constraints, proxy.prob_edge)
    return pdf(binomial_distribution, cost)
end


"""
N(c') from paper
"""
function N_cost_distribution(proxy::PaperProxy, cost::Integer)::Float64
    scale = 1 << proxy.num_qubits
    return P_cost_distribution(proxy, cost) * scale
end

"""
N(c'; d, c) from paper
"""
function N_cost_distance_distribution(proxy::PaperProxy,
        cost_1::Integer, distance::Integer, cost_2::Integer)::Float64

    sum = 0
    start_index = max(0, cost_1 + cost_2 - proxy.num_constraints)
    end_index = min(cost_1, cost_2)
    for common_constraints in start_index:end_index
        sum += prob_common_at_distance(proxy, common_constraints, cost_1, distance, cost_2)
    end

    p_cost = P_cost_distribution(proxy, cost_1)
    return (binomial(proxy.num_qubits, distance) / p_cost) * sum
end

"""
P(b, c'-b, c-b | d) from paper
"""
function prob_common_at_distance(proxy::PaperProxy, common_constraints::Integer,
        cost_1::Integer, distance::Integer, cost_2::Integer)::Float64

    #prob_same = (math.comb(proxy.num_constraints - distance, 2) + math.comb(distance, 2)) / math.comb(num_constraints, 2)
    prob_same = (binomial(proxy.num_qubits - distance, 2) + binomial(distance, 2)) / binomial(proxy.num_qubits, 2)
    prob_neither = prob_same / 2
    prob_both = prob_neither

    prob_one = (1 - prob_neither - prob_both) / 2
    probability_vec = SVector(prob_both, prob_one, prob_one, prob_neither)
    multinomial_distribution = Multinomial(proxy.num_constraints, probability_vec)
    k_vec = SVector(
        common_constraints,
        cost_1 - common_constraints,
        cost_2 - common_constraints,
        proxy.num_constraints + common_constraints - (cost_1 + cost_2)
    )
    return pdf(multinomial_distribution, k_vec)
end

