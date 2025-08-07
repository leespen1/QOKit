#%%  imports 
#%% imports 
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import qokit.maxcut as mc
from grips.QAOA_simulator import get_expectation, get_simulator, inverse_objective_function, QAOA_run
from grips.QAOA_proxy_interface import QAOA_proxy, QAOA_proxy_expectation
import grips.triangle_proxy as tpr
import grips.paper_proxy as ppr
import grips.normal_proxy as npr
import os  
import grips.real_distribution as rd 
from grips.triangle_proxy import TriangleProxy
from grips.QAOA_proxy_interface import QAOA_proxy_optimize_gamma_beta
from grips.scipy_additional_optimizers import spsa_for_scipy

#%%  Julia imports
from juliacall import Main as jl
jl.seval('''
using Pkg
Pkg.activate(joinpath(@__DIR__, "../julia"))
Pkg.instantiate()
using JuliaQAOA
''')

# %%
num_nodes = 3
edge_probability = 0.3
num_graphs = 2
graphs = [nx.erdos_renyi_graph(num_nodes, edge_probability) for _ in range(num_graphs)]

realdist = rd.get_homogeneous_distribution(graphs)/num_graphs
print(realdist)
realdist.shape



# %%
import numpy as np
from scipy.optimize import minimize

def mse_dist_loss(params, realdist, num_constraints, num_qubits):
    h_tweak_sub, hc_tweak_add, l_tweak_mul, r_tweak_mul = params
    
    proxy = TriangleProxy(
        num_constraints=num_constraints,
        num_qubits=num_qubits,
        h_tweak_sub=h_tweak_sub,
        hc_tweak_add=hc_tweak_add,
        l_tweak_mul=l_tweak_mul,
        r_tweak_mul=r_tweak_mul
    )
    
    predicted = np.zeros_like(realdist)
    
    # Loop over all cost_1, distance, cost_2
    for cost_1 in range(num_constraints+1):
        for distance in range(num_qubits+1):
            for cost_2 in range(num_constraints+1):
                predicted[cost_2, distance, cost_1] = proxy.N_cost_distance_distribution(cost_1, distance, cost_2)
    
    # Normalize both to sum to 1 (or same scale)
    predicted /= predicted.sum()
    realdist_norm = realdist / realdist.sum()
    
    # Compute MSE
    mse = np.mean((predicted - realdist_norm)**2)
    
    return mse

#%%
# Initial guess and bounds
initial_params = [0, 0, 1, 1]  
num_constraints = max((graphs[i].number_of_edges()) for i in range(num_graphs))  # Number of constraints -- +1 here caused error previously!
num_qubits = num_nodes  # Number of qubits
bounds = [
    (0, None),   # h_tweak_sub >= 0
    (-10, 10),   # hc_tweak_add can be small positive/negative
    (0.1, 20),   # l_tweak_mul > 0
    (0.1, 20),   # r_tweak_mul > 0
]

result = minimize(
    mse_dist_loss,
    initial_params,
    args=(realdist, num_constraints, num_qubits),
    method=spsa_for_scipy,  
    bounds=bounds,
    options={'maxiter': 1}
)

fitted_params = result.x
print("Fitted parameters:", fitted_params)
fitted_proxy = TriangleProxy(
    num_constraints=num_constraints,
    num_qubits=num_qubits,
    h_tweak_sub=fitted_params[0],
    hc_tweak_add=fitted_params[1],
    l_tweak_mul=fitted_params[2],
    r_tweak_mul=fitted_params[3]
)


#%%
gammas = np.linspace(0, np.pi, 10)  # Gamma values for QAOA
betas = np.linspace(0, np.pi, 10)  # Beta values for QAOA
fitted_triangle_results = QAOA_proxy(fitted_proxy, gammas, betas)
print("Fitted Proxy Results:", fitted_triangle_results)

final_amplitudes = fitted_triangle_results[-1]
expectation = QAOA_proxy_expectation(fitted_proxy, final_amplitudes)
print("Expectation value:", expectation)

# %%
# Compare MSE for initial and fitted parameters
initial_mse = mse_dist_loss(initial_params, realdist, num_constraints, num_qubits)
fitted_mse = mse_dist_loss(fitted_params, realdist, num_constraints, num_qubits)
print(f"Initial MSE: {initial_mse}")
print(f"Fitted MSE: {fitted_mse}")

# %%
# Compute results with initial parameters
initial_proxy = TriangleProxy(
    num_constraints=num_constraints,
    num_qubits=num_qubits,
    h_tweak_sub=initial_params[0],
    hc_tweak_add=initial_params[1],
    l_tweak_mul=initial_params[2],
    r_tweak_mul=initial_params[3]
)

initial_triangle_results = QAOA_proxy(initial_proxy, gammas, betas)
print("Initial Proxy Results:", initial_triangle_results)

# initial_final_amplitudes = initial_triangle_results[-1]
# initial_expectation = QAOA_proxy_expectation(initial_proxy, initial_final_amplitudes)
# print("Initial Expectation value:", initial_expectation)

# # Compare expectations
# print(f"Initial Expectation: {initial_expectation}")
# print(f"Fitted Expectation: {expectation}")


# %% doing this correctly now, I think: 
'''
-defining initial gammas and betas
-finding best gamma and beta according to initial versus fitted proxy
-using these to run QAOA and comparing expectations
'''

gamma_0 = np.array([0.1])
beta_0 = np.array([0.1])
init_result = QAOA_proxy_optimize_gamma_beta(initial_proxy, gamma_0, beta_0, optimizer_options={'maxiter': 1})
gamma_init = init_result["gamma"]
beta_init = init_result["beta"]



fitted_result = QAOA_proxy_optimize_gamma_beta(fitted_proxy, gamma_0, beta_0)
gamma_fitted = fitted_result["gamma"]
beta_fitted = fitted_result["beta"]

num_graphs = 2
graphs = [nx.erdos_renyi_graph(num_nodes, edge_probability) for _ in range(num_graphs)]


#get QAOA expectations of inverse objective function for initial and fitted proxies
initial_expectations = []
fitted_expectations = []
for graph in graphs:
    ising_model = mc.get_maxcut_terms(graph)
    N = graph.number_of_nodes()
    sim = get_simulator(N, ising_model)
    p = 1
    mixer = "x"
    expectations = []
    overlaps = []
    initres = QAOA_run(
        ising_model=ising_model,
        N=N,
        p=p,
        init_gamma=gamma_init,
        init_beta=beta_init,
        optimizer_method="COBYLA",
        optimizer_options=None,
        mixer=mixer,
        expectations=expectations,
        overlaps=overlaps
    )
    print("Initial Proxy Result:", initres)

    # For fitted parameters
    expectations = []
    overlaps = []
    fitres = QAOA_run(
        ising_model=ising_model,
        N=N,
        p=p,
        init_gamma=gamma_fitted,
        init_beta=beta_fitted,
        optimizer_method="COBYLA",
        optimizer_options=None,
        mixer=mixer,
        expectations=expectations,
        overlaps=overlaps
    )
    print("Fitted Proxy Result:", fitres)

    # Store the expectations
    initial_expectations.append(initres["expectation"])
    fitted_expectations.append(fitres["expectation"])

# Print the mean expectations
print("Initial Proxy Mean Expectation:", np.mean(initial_expectations))
print("Fitted Proxy Mean Expectation:", np.mean(fitted_expectations))


# %%
