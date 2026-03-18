"""
This file implements the interface for approximating QAOA using a homogenous
proxy. The main functions are:
- QAOA_proxy: *runs* the homogeneous proxy of the QAOA circuit.
- QAOA_proxy_expectation: computes the *expectation value* of a proxy state.
- QAOA_proxy_optimize_gamma_beta: optimizes the parameters gamma and beta of the
                                  QAOA circuit using the homogeneous proxy.

Several function also have an @njit version. That is done because the proxies
that call external libraries cannot be compiled using Numba, so we can't always
use @njit. To deal with this, we basically have to write the function twice,
once without @njit and once with @njit, then have the version without @njit
call the @njit version if the proxy can be compiled.
"""
from numba import njit, jit
import numpy as np
import typing, time, scipy, os
from scipy.stats import binom, multinomial


# Make the functions written in Julia available (call with `jl.function_name`)
# (by having this in __init__.py, should be able to simply do from grips import jl)
from juliacall import Main as jl
grips_dir = os.path.dirname(os.path.abspath(__file__))
julia_project_dir = os.path.normpath(os.path.join(grips_dir, "../julia"))
jl.seval(f'''
using Pkg
Pkg.activate("{julia_project_dir}")
try
    using JuliaQAOA
catch e
    println("Encountered error during 'using JuliaQAOA', instantiating environment at {julia_project_dir}")
    Pkg.instantiate()
    using JuliaQAOA
end
#try
#    using JuliaQAOA
#catch e
#    if isa(e, LoadError)
#        println("Encountered error during 'using JuliaQAOA', instantiating environment at {julia_project_dir}")
#        Pkg.instantiate()
#        using JuliaQAOA
#    else
#        rethrow(e)
#    end
#end
''')




def QAOA_proxy(
    proxy,
    gammas: np.ndarray, # 1D array
    betas: np.ndarray,  # 2D array
) -> np.ndarray:
    """
    Run the homogeneous proxy of the QAOA circuit with parameters gammas (cost) 
    and betas (mixer).

    Return a complex matrix where each row is the set of proxy amplitudes
    after each "layer" of the QAOA circuit (i.e., the homogeneous proxy of the 
    state vector).
    """

    # If proxy is a NJIT-compiled obj, use njit version
    if type(proxy).__module__ == 'numba.experimental.jitclass.boxing':
        amplitude_proxies = QAOA_proxy_njit(proxy, gammas, betas)
        return amplitude_proxies
    # If proxy is is a Juliacall obj, use julia version
    elif type(proxy).__module__ == 'juliacall':
        # Convert result to numpy array, and take transpose (Because Julia is column-major, while Numpy is row-major)
        amplitude_proxies = jl.QAOA_proxy(proxy, gammas, betas).to_numpy().T
        return amplitude_proxies

    # Make sure gammas and betas are vectors, which have same shape
    assert gammas.shape == betas.shape
    assert len(gammas.shape) == 1
    num_QAOA_layers = len(gammas)

    num_unique_costs = proxy.num_constraints + 1
    # The p-th row gives the probability amplitudes for each cost after the p-th
    # layer of the QAOA homogenous proxy
    amplitude_proxies = np.zeros(
        (num_QAOA_layers + 1, num_unique_costs),
        dtype=np.complex128
    ) # (p+1, num_costs) needs to be a tuple, not a list, in order to play nicely with numba. Also, dtype must be made more concrete (complex128 instead of complex)

    # Set up initial amplitudes
    init_amplitude = np.sqrt(1 / (1 << proxy.num_qubits))
    for i in range(num_unique_costs):
        amplitude_proxies[0][i] = init_amplitude

    # Run the QAOA proxy, get probability amplitudes after each layer
    for current_depth in range(1, num_QAOA_layers + 1):
        for cost_1 in range(num_unique_costs):
            prev_amplitude_proxies = amplitude_proxies[current_depth - 1]
            amplitude_proxies[current_depth][cost_1] = compute_amplitude_sum(
                proxy,
                prev_amplitude_proxies,
                gammas[current_depth - 1],
                betas[current_depth - 1],
                cost_1,
            )

    return amplitude_proxies


def QAOA_proxy_expectation(
    proxy,
    amplitude_proxies: np.ndarray, # 1D array
    drop_first_N_costs: int = 0 # Optionally, ignore the lowest costs 
) -> float:
    """
    Given a set of proxy amplitudes (i.e. the homomgeneous proxy of the state
    vector.), compute the expectation value associated with those amplitudes.
    """

    # If proxy is a NJIT-compiled obj, use njit version
    if type(proxy).__module__ == 'numba.experimental.jitclass.boxing':
        proxy_expectation_value = QAOA_proxy_expectation_njit(
            proxy, amplitude_proxies, drop_first_N_costs
        )
        return proxy_expectation_value
    # If proxy is is a Juliacall obj, use julia version
    elif type(proxy).__module__ == 'juliacall':
        # Convert result to numpy array, and take transpose (Because Julia is column-major, while Numpy is row-major)
        proxy_expectation_value = jl.QAOA_proxy_expectation(
            proxy, amplitude_proxies, drop_first_N_costs
        )
        return proxy_expectation_value
            
    assert len(amplitude_proxies.shape) == 1, f"Bad array size {amplitude_proxies.shape}."

    proxy_expectation_value = 0
    num_unique_costs = amplitude_proxies.shape[0]

    for cost in range(drop_first_N_costs, num_unique_costs):
        proxy_expectation_value += cost * proxy.N_cost_distribution(cost) * (abs(amplitude_proxies[cost]) ** 2)

    return proxy_expectation_value


def QAOA_proxy_expectation_from_gamma_beta(
    proxy,
    gammas: np.ndarray, # 1D array
    betas: np.ndarray,  # 2D array
) -> float:
    """
    Convenience function that runs the QAOA proxy and then computes the
    expectation value, given parameters gammas and betas.
    """

    amplitude_proxies = QAOA_proxy(
        proxy,
        gammas,
        betas,
    )
    final_amplitude_proxies = amplitude_proxies[-1]

    expectation = QAOA_proxy_expectation(
        proxy, final_amplitude_proxies
    )

    return expectation



def QAOA_proxy_optimize_gamma_beta(
    proxy,
    init_gamma: np.ndarray,
    init_beta: np.ndarray,
    optimizer_method: str = "COBYLA",
    optimizer_options: dict | None = None,
    expectations: list[np.ndarray] | None = None,
) -> dict:
    """
    Using a homogeneous proxy of the QAOA circuit, optimize gamma and beta to
    find values which maximize the expectation value of the final state of the
    circuit.

    Hopefully, these parameters will also result in the real QAOA circuit
    resulting in a state with high expectation value.
    """

    assert init_gamma.shape == init_beta.shape
    assert len(init_gamma.shape) == 1
    num_QAOA_layers = len(init_gamma)

    init_freq = np.hstack([init_gamma, init_beta])

    start_time = time.time()
    result = scipy.optimize.minimize(
        inverse_proxy_objective_function(
            proxy,
            num_QAOA_layers,
            expectations
        ),
        init_freq,
        args=(),
        method=optimizer_method,
        options=optimizer_options,
    )
    # the above returns a scipy optimization result object that has multiple attributes
    # result.x gives the optimal solutionsol.success #bool whether algorithm succeeded
    # result.message #message of why algorithms terminated
    # result.nfev is number of iterations used (here, number of QAOA calls)
    end_time = time.time()

    def make_time_relative(input: tuple[float, float]) -> tuple[float, float]:
        time, x = input
        return (time - start_time, x)

    if expectations is not None:
        expectations = list(map(make_time_relative, expectations))

    # Run a final time, to get the expectation with the optimized parameters
    gamma, beta = result.x[:num_QAOA_layers], result.x[num_QAOA_layers:]
    expectation = -result.fun # TODO check that this is correct

    return {
        "gamma": gamma,
        "beta": beta,
        "expectation": expectation,
        "runtime": end_time - start_time,  # measured in seconds
        "num_QAOA_calls": result.nfev,  # Calls to the proxy
        "classical_opt_success": result.success,
        "scipy_opt_message": result.message,
    }



def compute_amplitude_sum(
    proxy,
    prev_amplitudes: np.ndarray, # Q_l(c) for each unique 
    gamma: float, # QAOA parameter gamma, how long to apply cost hamiltonian
    beta: float,  # QAOA parameter beta, how long to apply mixer hamiltonian
    cost_1: int,  # c', the cost whose probability amplitude we want to compute
) -> complex:
    """
    Given parameters gamma and beta, and the set of probability amplitudes
    Q_l(c) for each unique cost c after performing l "layers" of the QAOA
    homogenous proxy, and a particular cost c', return Q_l+1(c'), which is the
    amplitude of cost c' after the l+1-th layer of the QAOA homogeneous proxy.

    See the for-loop of Algorithm 1 in parameter-setting paper.
    """

    assert len(prev_amplitudes.shape) == 1, f"Bad array size {prev_amplitudes.shape}."

    # If proxy is a NJIT-compiled obj, use njit version
    if type(proxy).__module__ == 'numba.experimental.jitclass.boxing':
        amp_sum = compute_amplitude_sum_njit(
            proxy, prev_amplitudes, gamma, beta, cost_1
        )
        return amp_sum

    # Non-compiled version
    amp_sum = 0
    num_unique_costs = prev_amplitudes.shape[0]

    for cost_2 in range(num_unique_costs):
        for distance in range(proxy.num_qubits + 1):
            beta_factor = (np.cos(beta) ** (proxy.num_qubits - distance)) * ((-1j * np.sin(beta)) ** distance)
            gamma_factor = np.exp(-0.5j * gamma * cost_2)
            num_costs_at_distance = proxy.N_cost_distance_distribution(
                cost_1, distance, cost_2
            )
            amp_sum += beta_factor * gamma_factor * prev_amplitudes[cost_2] * num_costs_at_distance

    return amp_sum



@njit
def compute_amplitude_sum_njit(
    proxy,
    prev_amplitudes: np.ndarray, # Q_l(c) for each unique
    gamma: float, # QAOA parameter gamma, how long to apply cost hamiltonian
    beta: float,  # QAOA parameter beta, how long to apply mixer hamiltonian
    cost_1: int,  # c', the cost whose probability amplitude we want to compute
) -> complex:
    """
    Given parameters gamma and beta, and the set of probability amplitudes
    Q_l(c) for each unique cost c after performing l "layers" of the QAOA
    homogenous proxy, and a particular cost c', return Q_l+1(c'), which is the
    amplitude of cost c' after the l+1-th layer of the QAOA homogeneous proxy.

    Convention: uses exp(-iγc/2), matching QOKit's `exp(-0.5j * gamma * hc_diag)`.

    TODO: This function is not NJIT-compiled because the method
    proxy.N_cost_distribution is not necessarily NJIT-compiled.
    """

    assert len(prev_amplitudes.shape) == 1, f"Bad array size {prev_amplitudes.shape}."

    amp_sum = 0
    num_unique_costs = prev_amplitudes.shape[0]

    for cost_2 in range(num_unique_costs):
        for distance in range(proxy.num_qubits + 1):
            beta_factor = (np.cos(beta) ** (proxy.num_qubits - distance)) * ((-1j * np.sin(beta)) ** distance)
            gamma_factor = np.exp(-0.5j * gamma * cost_2)
            num_costs_at_distance = proxy.N_cost_distance_distribution(
                cost_1, distance, cost_2
            )
            amp_sum += beta_factor * gamma_factor * prev_amplitudes[cost_2] * num_costs_at_distance

    return amp_sum



@njit
def QAOA_proxy_njit(
    proxy,
    gammas: np.ndarray, # 1D array
    betas: np.ndarray,  # 2D array
) -> np.ndarray:

    # Make sure gammas and betas are vectors, which have same shape
    assert gammas.shape == betas.shape
    assert len(gammas.shape) == 1
    num_QAOA_layers = len(gammas)

    num_unique_costs = proxy.num_constraints + 1
    # The p-th row gives the probability amplitudes for each cost after the p-th
    # layer of the QAOA homogenous proxy
    amplitude_proxies = np.zeros(
        (num_QAOA_layers + 1, num_unique_costs),
        dtype=np.complex128
    ) # (p+1, num_costs) needs to be a tuple, not a list, in order to play nicely with numba. Also, dtype must be made more concrete (complex128 instead of complex)

    # Set up initial amplitudes
    init_amplitude = np.sqrt(1 / (1 << proxy.num_qubits))
    for i in range(num_unique_costs):
        amplitude_proxies[0][i] = init_amplitude

    # Run the QAOA proxy, get probability amplitudes after each layer
    for current_depth in range(1, num_QAOA_layers + 1):
        for cost_1 in range(num_unique_costs):
            prev_amplitude_proxies = amplitude_proxies[current_depth - 1]
            amplitude_proxies[current_depth][cost_1] = compute_amplitude_sum_njit(
                proxy,
                prev_amplitude_proxies,
                gammas[current_depth - 1],
                betas[current_depth - 1],
                cost_1,
            )

    return amplitude_proxies



@njit
def QAOA_proxy_expectation_njit(
    proxy,
    amplitude_proxies: np.ndarray,
    drop_first_N_costs: int = 0 # Optionally, ignore the lowest costs 
) -> float:

    assert len(amplitude_proxies.shape) == 1, f"Bad array size {amplitude_proxies.shape}."

    proxy_expectation_value = 0
    num_unique_costs = amplitude_proxies.shape[0]

    for cost in range(drop_first_N_costs, num_unique_costs):
        proxy_expectation_value += cost * proxy.N_cost_distribution(cost) * (abs(amplitude_proxies[cost]) ** 2)

    return proxy_expectation_value



def inverse_proxy_objective_function(
    proxy,
    num_QAOA_layers: int,
    expectations: list[np.ndarray] | None
) -> typing.Callable:
    """
    Sets ups the objective function to be optimized in
    `QAOA_proxy_optimize_gamma_beta`.

    Inverse here is the additive inverse, not multiplicative inverse.

    I.e. we are taking -x, not 1/x.

    Scipy allows us to minimize functions, but we want to maximize our cost
    (the number of edges cut), so we have to flip the +/- sign to turn our
    problem into a minimization problem.
    """

    # If proxy is a NJIT-compiled obj, use njit version
    if type(proxy).__module__ == 'numba.experimental.jitclass.boxing':
        def inverse_objective(*args) -> float:
            gamma, beta = args[0][:num_QAOA_layers], args[0][num_QAOA_layers:]
             
            amplitude_proxies = QAOA_proxy_njit(proxy, gamma, beta)
            final_amplitude_proxies = amplitude_proxies[-1]
            expectation = QAOA_proxy_expectation_njit(proxy, final_amplitude_proxies)

            current_time = time.time()

            if expectations is not None:
                expectations.append((current_time, expectation))

            return -expectation

        return inverse_objective
    elif type(proxy).__module__ == 'juliacall':
        inverse_objective = jl.inverse_proxy_objective_function(
            proxy, num_QAOA_layers, expectations
        )
        return inverse_objective

    def inverse_objective(*args) -> float:
        gamma, beta = args[0][:num_QAOA_layers], args[0][num_QAOA_layers:]
         
        amplitude_proxies = QAOA_proxy(
            proxy,
            gamma,
            beta,
        )
        final_amplitude_proxies = amplitude_proxies[-1]

        expectation = QAOA_proxy_expectation(
            proxy,
            final_amplitude_proxies
        )
        current_time = time.time()

        if expectations is not None:
            expectations.append((current_time, expectation))

        return -expectation

    return inverse_objective





import numpy as np
import networkx as nx
from scipy.stats import norm, beta, skew, kurtosis
from collections import defaultdict
import math

class CostDistribution:
    """
    A callable wrapper for the estimated distribution P(c').
    Now supports advanced continuous approximations.
    """
    def __init__(self, method, probability_map=None, params=None, dist_type=None):
        self.method = method
        self.probability_map = probability_map # For Wang-Landau
        self.params = params # For Moment-Matching
        self.dist_type = dist_type # 'gaussian', 'beta', 'edgeworth'
        
        if probability_map:
            self.support = set(probability_map.keys())
        
    def __call__(self, c):
        if self.method == 'wang_landau':
            # Discrete lookup
            if isinstance(c, (int, np.integer)):
                return self.probability_map.get(c, 0.0)
            c_round = int(round(c))
            if abs(c - c_round) < 1e-5:
                return self.probability_map.get(c_round, 0.0)
            return 0.0
            
        elif self.method == 'moment_matching':
            if self.dist_type == 'gaussian':
                mu, sigma = self.params
                return norm.pdf(c, loc=mu, scale=sigma)
            
            elif self.dist_type == 'beta':
                # params = (a, b, loc, scale)
                return beta.pdf(c, *self.params)
            
            elif self.dist_type == 'edgeworth':
                # Edgeworth expansion A (Gram-Charlier) around Gaussian
                # params = (mu, sigma, skewness, kurtosis)
                mu, sigma, S, K = self.params
                z = (c - mu) / sigma
                phi = norm.pdf(z)
                # Expansion up to 4th moment (Kurtosis)
                # H3(z) = z^3 - 3z
                # H4(z) = z^4 - 6z^2 + 3
                H3 = z**3 - 3*z
                H4 = z**4 - 6*z**2 + 3
                
                term1 = (S / 6) * H3
                term2 = (K / 24) * H4
                term3 = (S**2 / 72) * (z**6 - 15*z**4 + 45*z**2 - 15) # higher order skew correction
                
                # p(z) = phi(z) * (1 + ...)
                # Note: We omit term3 for stability if desired, but standard Edgeworth includes it.
                # We return density in x space: 1/sigma * p(z)
                density = (1 / sigma) * phi * (1 + term1 + term2)
                return max(0.0, density) # Ensure non-negative

def _get_cut_size_vectorized(adj_matrix, bitstrings):
    """Computes cut sizes for a batch of bitstrings efficiently."""
    spins = 2 * bitstrings - 1
    interactions = (spins @ adj_matrix) * spins 
    interaction_sum = interactions.sum(axis=1)
    num_edges = adj_matrix.sum() / 2
    return 0.5 * (num_edges - 0.5 * interaction_sum)

def _run_wang_landau(graph, steps):
    """Runs Wang-Landau sampling."""
    n = graph.number_of_nodes()
    adj = nx.to_dict_of_lists(graph)
    x = np.random.randint(0, 2, n)
    
    current_cost = 0
    for u, v in graph.edges():
        if x[u] != x[v]: current_cost += 1
            
    log_g = defaultdict(float)
    log_f = 1.0
    decay_interval = max(100, steps // 10)
    
    for t in range(steps):
        flip_node = np.random.randint(0, n)
        delta = 0
        original_val = x[flip_node]
        for neighbor in adj[flip_node]:
            if original_val != x[neighbor]: delta -= 1
            else: delta += 1
        
        candidate_cost = current_cost + delta
        if log_g[current_cost] - log_g[candidate_cost] >= np.log(np.random.rand()):
            current_cost = candidate_cost
            x[flip_node] = 1 - x[flip_node]
        
        log_g[current_cost] += log_f
        if (t + 1) % decay_interval == 0: log_f /= 2.0
            
    return log_g

def estimate_P_for_graphs(graphs, method='moment_matching', num_samples=10000, dist_type='edgeworth'):
    """
    Estimates P(c') with improved moment matching options.
    
    Args:
        method (str): 'wang_landau' or 'moment_matching'
        dist_type (str): Used if method='moment_matching'. 
                         Options: 'gaussian', 'beta', 'edgeworth'.
                         'beta' is recommended for bounded graph problems.
    """
    if not graphs: raise ValueError("Graph list is empty.")
    num_graphs = len(graphs)
    samples_per_graph = max(1, num_samples // num_graphs)
    
    if method == 'moment_matching':
        all_costs = []
        for G in graphs:
            n = G.number_of_nodes()
            adj = nx.to_scipy_sparse_array(G)
            bitstrings = np.random.randint(0, 2, size=(samples_per_graph, n))
            all_costs.extend(_get_cut_size_vectorized(adj, bitstrings))
        
        all_costs = np.array(all_costs)
        
        if dist_type == 'gaussian':
            mu, sigma = norm.fit(all_costs)
            return CostDistribution('moment_matching', params=(mu, sigma), dist_type='gaussian')
            
        elif dist_type == 'beta':
            # Beta fits 4 parameters: alpha, beta, loc (min), scale (range)
            # It handles skewness and finite support very well.
            a, b, loc, scale = beta.fit(all_costs)
            return CostDistribution('moment_matching', params=(a, b, loc, scale), dist_type='beta')
            
        elif dist_type == 'edgeworth':
            mu = np.mean(all_costs)
            sigma = np.std(all_costs)
            S = skew(all_costs)
            K = kurtosis(all_costs) # Fisher kurtosis (normal=0)
            return CostDistribution('moment_matching', params=(mu, sigma, S, K), dist_type='edgeworth')
            
    elif method == 'wang_landau':
        # Wang-Landau logic (same as before, aggregating histograms)
        combined_prob_map = defaultdict(float)
        for G in graphs:
            log_g = _run_wang_landau(G, samples_per_graph)
            max_log = max(log_g.values())
            probs = {c: np.exp(val - max_log) for c, val in log_g.items()}
            norm_factor = sum(probs.values())
            for c, val in probs.items():
                combined_prob_map[c] += (val / norm_factor) / num_graphs
        return CostDistribution('wang_landau', probability_map=dict(combined_prob_map))

    else:
        raise ValueError(f"Unknown method: {method}")

