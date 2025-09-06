#%% imports
print("Importing python packages/functions ...")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import qokit.maxcut as mc
import os
import grips
from grips import (
    get_simulator, QAOA_run, QAOA_proxy, QAOA_proxy_expectation, QAOA_proxy_optimize_gamma_beta,
    get_homogeneous_distribution,
    TriangleProxy, NormalProxy, PaperProxy,
    inverse_objective_function, get_expectation,  
    maxcut, maxcut_approx_ratio, spsa_for_scipy,
    plot_distribution_lines_all
)
from grips.sendai_opt import fit_proxy_to_real
from scipy.optimize import minimize
from scipy.optimize import dual_annealing
print("Finished importing python packages/functions!")

#%%  Julia imports
print("Importing Julia functions ...")
from juliacall import Main as jl
jl.seval('''
using Pkg
Pkg.activate(joinpath(@__DIR__, "../julia"))
Pkg.instantiate()
using JuliaQAOA
''')
print("Finished importing Julia functions!")

# %% Set up a single graph, get its statistical homogeneous distribution
# 3 nodes, edge probability 0.3, seed=4 results in a graph with 2 edges
print("\nSetting up a single graph, getting its statistical homogeneous distribution ...")
num_nodes = 3
edge_probability = 0.3
graph = nx.erdos_renyi_graph(num_nodes, edge_probability)

#I think graphs with no edges were causing Nans! -PK
if graph.number_of_edges() == 0 and num_nodes > 1:
    nodes = list(graph.nodes())
    u, v = np.random.choice(nodes, 2, replace=False)
    graph.add_edge(u,v)
    print("Graph had no edges, added a random edge.")

realdist = rd.get_homogeneous_distribution(graph)
realdist = np.nan_to_num(realdist) #pad with zeros instead in case of nans
print(realdist)

realdist.shape

print("The graph:")
nx.draw(graph, with_labels=True, node_color="skyblue", node_size=2000, font_size=14, font_weight="bold")
print("The homogeneous distribution:")
plot_distribution_lines_all(homodist, "Averaged Homogeneous Distribution for Random Graph")

print("Finished setting up graph and getting homogeneous distribution!")
#I think graphs with no edges were causing Nans! -PK
if graph.number_of_edges() == 0 and num_nodes > 1:
    nodes = list(graph.nodes())
    u, v = np.random.choice(nodes, 2, replace=False)
    graph.add_edge(u,v)
    print("Graph had no edges, added a random edge.")

realdist = get_homogeneous_distribution(graph)
realdist = np.nan_to_num(realdist) #pad with zeros instead in case of nans
print(realdist)

realdist.shape


#%% Show initial proxy distribution
initial_params = [100, 0, 1, 1]  
num_constraints = graph.number_of_edges() # Number of constraints -- +1 here caused error previously!
num_qubits = num_nodes  # Number of qubits
proxy = TriangleProxy(num_constraints, num_qubits, *initial_params)
initial_triangle_homodist = grips.get_homogeneous_distribution_from_proxy(proxy)
plot_distribution_lines_all(initial_triangle_homodist, f"Initial Triangle Proxy Distribution {initial_params}")

# %% Define mean-squared error loss function for fitting triangle proxy to statistical homogeneous distribution
print("\nDefining triangle proxy loss function ... ")

# def mse_dist_loss(params, realdist, num_constraints, num_qubits):
#     h_tweak_sub, hc_tweak_add, l_tweak_mul, r_tweak_mul = params
    
#     proxy = TriangleProxy(
#         num_constraints=num_constraints,
#         num_qubits=num_qubits,
#         h_tweak_sub=h_tweak_sub,
#         hc_tweak_add=hc_tweak_add,
#         l_tweak_mul=l_tweak_mul,
#         r_tweak_mul=r_tweak_mul
#     )
    
#     predicted = np.zeros_like(realdist)
    
#     # Loop over all cost_1, distance, cost_2
#     #note: this is parallelizable
#     for cost_1 in range(num_constraints+1):
#         for distance in range(num_qubits+1):
#             for cost_2 in range(num_constraints+1):
#                 predicted[cost_2, distance, cost_1] = proxy.N_cost_distance_distribution(cost_1, distance, cost_2)
    
#     # Normalize both to sum to 1 (or same scale)
#     predicted /= predicted.sum()
#     realdist_norm = realdist / realdist.sum()
    
#     # Compute MSE
#     mse = np.mean((predicted - realdist_norm)**2)
    
#     return mse

#%%
print("\nFittting triangle proxy to real distribution ...")
# Initial guess and bounds
use_small_bounds = True
initial_params = [0, 0, 1, 1]  
num_constraints = graph.number_of_edges() # Number of constraints -- +1 here caused error previously!
num_qubits = num_nodes  # Number of qubits
bounds = [
    (0, None),   # h_tweak_sub >= 0
    (-10, 10),   # hc_tweak_add can be small positive/negative
    (0.1, 20),   # l_tweak_mul > 0
    (0.1, 20),   # r_tweak_mul > 0
]

# result = minimize(
#     mse_dist_loss,
#     initial_params,
#     args=(homodist, num_constraints, num_qubits),
#     method='Nelder-Mead',  
#     bounds=bounds,
#     options={'maxiter': 10000, 'epsilon': 0.0001}
# )

bounds[0] = (0,10) #dual annealing needs non-None bounds 

small_bounds = [
    (0, 3),   # h_tweak_sub >= 0
    (-3, 3),   # hc_tweak_add can be small positive/negative
    (0.3, 5),   # l_tweak_mul > 0
    (0.3, 5),   # r_tweak_mul > 0
]
if use_small_bounds:
    bounds = small_bounds

# result = dual_annealing(
#     mse_dist_loss,
#     bounds=bounds,
#     args=(realdist, num_constraints, num_qubits),
#     maxiter=50000  #this is almost definitely too many its but works for now :) 
# )
# result = dual_annealing(
#     mse_dist_loss,
#     bounds=bounds,
#     args=(realdist, num_constraints, num_qubits),
#     maxiter=50000  #this is almost definitely too many its but works for now :) 
# )

# fitted_params = result.x

startproxy = TriangleProxy(
    num_constraints=num_constraints,
    num_qubits=num_qubits,
    h_tweak_sub=initial_params[0],
    hc_tweak_add=initial_params[1],
    l_tweak_mul=initial_params[2],
    r_tweak_mul=initial_params[3]
)
fitted_params, _ = fit_proxy_to_real(startproxy, realdist, initial_params, bounds, num_constraints,\
                      num_qubits, max_iter = 1000)

# fitted_params = result.x
initial_proxy = TriangleProxy(
    num_constraints=num_constraints,
    num_qubits=num_qubits,
    h_tweak_sub=initial_params[0],
    hc_tweak_add=initial_params[1],
    l_tweak_mul=initial_params[2],
    r_tweak_mul=initial_params[3]
)

#the fit proxy to real modified the proxy in-place, so starting a new instance.
startproxy_for_opt = TriangleProxy(
    num_constraints=num_constraints,
    num_qubits=num_qubits,
    h_tweak_sub=initial_params[0],
    hc_tweak_add=initial_params[1],
    l_tweak_mul=initial_params[2],
    r_tweak_mul=initial_params[3]
)
fitted_params, _ = fit_proxy_to_real(startproxy_for_opt, realdist, initial_params, bounds, num_constraints,\
                      num_qubits, max_iter = 1000)

print("Fitted parameters:", fitted_params)
fitted_proxy = TriangleProxy(
    num_constraints=num_constraints,
    num_qubits=num_qubits,
    h_tweak_sub=fitted_params[0],
    hc_tweak_add=fitted_params[1],
    l_tweak_mul=fitted_params[2],
    r_tweak_mul=fitted_params[3]
)
fitted_proxy_homodist = grips.get_homogeneous_distribution_from_proxy(fitted_proxy)
plot_distribution_lines_all(fitted_proxy_homodist, f"Fitted Proxy Homogeneous Distribution for Random Graph {fitted_params}")

print("Finished fitting triangle proxy to homogeneous distribution!")

# %% Compare MSE for initial and fitted parameters
print("\nComparing MSE loss function for initial and fitted parameters ...")
# Compare MSE for initial and fitted parameters
initial_mse = grips.distribution_mean_squared_error(initial_proxy, realdist)
fitted_mse = grips.distribution_mean_squared_error(fitted_proxy, realdist)
print(f"Initial MSE: {initial_mse}")
print(f"Fitted MSE: {fitted_mse}")
print("Finished comparing MSE loss function for initial and fitted parameters!")


#%% Run QAOA with fitted proxy
print("\nRunning QAOA with fitted proxy (but not tuned gamma/beta)...")
gammas = np.linspace(0, np.pi, 10)  # Gamma values for QAOA
betas = np.linspace(0, np.pi, 10)  # Beta values for QAOA
fitted_triangle_results = QAOA_proxy(fitted_proxy, gammas, betas)
print("Fitted Proxy Results:", fitted_triangle_results)

final_amplitudes = fitted_triangle_results[-1]
expectation = QAOA_proxy_expectation(fitted_proxy, final_amplitudes)
print("Expectation value:", expectation)


# %% Run QAOA with initial, unfitted proxy
print("\nRunning QAOA with initial, unfitted proxy ...")

# %%
# Compute results with initial parameters
# initial_proxy = TriangleProxy(
#     num_constraints=num_constraints,
#     num_qubits=num_qubits,
#     h_tweak_sub=initial_params[0],
#     hc_tweak_add=initial_params[1],
#     l_tweak_mul=initial_params[2],
#     r_tweak_mul=initial_params[3]
# )

# # Compare MSE for initial and fitted parameters
# initial_mse = mse_dist_loss(initial_proxy, realdist, num_constraints)
# fitted_mse = mse_dist_loss(fitted_proxy, realdist, num_constraints)
# print(f"Initial MSE: {initial_mse}")
# print(f"Fitted MSE: {fitted_mse}")


# Compare MSE for initial and fitted parameters
initial_mse = mse_dist_loss(initial_proxy, realdist, num_constraints, num_qubits)
fitted_mse = mse_dist_loss(fitted_proxy, realdist, num_constraints, num_qubits)
print(f"Initial MSE: {initial_mse}")
print(f"Fitted MSE: {fitted_mse}")


initial_triangle_results = QAOA_proxy(initial_proxy, gammas, betas)
print("Initial Proxy Results:", initial_triangle_results)

# initial_final_amplitudes = initial_triangle_results[-1]
# initial_expectation = QAOA_proxy_expectation(initial_proxy, initial_final_amplitudes)
# print("Initial Expectation value:", initial_expectation)

# # Compare expectations
# print(f"Initial Expectation: {initial_expectation}")
# print(f"Fitted Expectation: {expectation}")
print("Finished running QAOA with initial, unfitted proxy!")


# Comparing QAOA results for initial versus fitted proxy
'''
-defining initial gammas and betas
-finding best gamma and beta according to initial versus fitted proxy
-using these to run QAOA and comparing expectations
'''

gamma_0 = np.array([0.1])
beta_0 = np.array([0.1])
init_result = QAOA_proxy_optimize_gamma_beta(initial_proxy, gamma_0, beta_0, optimizer_method = 'Nelder-Mead', optimizer_options={'maxiter': 1})
gamma_init = init_result["gamma"]
beta_init = init_result["beta"]



fitted_result = QAOA_proxy_optimize_gamma_beta(fitted_proxy, gamma_0, beta_0, optimizer_method = 'Nelder-Mead', 
                                               optimizer_options = {'maxiter': 10000, 'epsilon': 0.00001})
gamma_fitted = fitted_result["gamma"]
beta_fitted = fitted_result["beta"]

graph = nx.erdos_renyi_graph(num_nodes, edge_probability)


#%% get QAOA expectations of inverse objective function for initial and fitted proxies
initial_expectations = []
fitted_expectations = []

initial_approx_ratios = []
fitted_approx_ratios = []

ising_model = mc.get_maxcut_terms(graph)
N = graph.number_of_nodes()
sim = get_simulator(N, ising_model)
p = 1
mixer = "x"
expectations = []
overlaps = []

'''
we want the QAOA expectations of the fitted gammas and bets from init proxy versus fitted proxy, 
not using them as starter values to optimize. What is currently commented out is doing the latter, 
which we don't want, but may want to look at later too. 
'''

# initres = QAOA_run(
#     ising_model=ising_model,
#     N=N,
#     p=p,
#     init_gamma=gamma_init,
#     init_beta=beta_init,
#     optimizer_method="COBYLA",
#     optimizer_options=None,
#     mixer=mixer,
#     expectations=expectations,
#     overlaps=overlaps
# )
inv_obj = inverse_objective_function(ising_model, N, p, mixer, None, None)
initres = inv_obj(np.hstack([gamma_init, beta_init]))
print("Initial Proxy Result (negative expectation):", initres)

# For fitted parameters
expectations = []
overlaps = []
# fitres = QAOA_run(
#     ising_model=ising_model,
#     N=N,
#     p=p,
#     init_gamma=gamma_fitted,
#     init_beta=beta_fitted,
#     optimizer_method='COBYLA',
#     optimizer_options=None,
#     mixer=mixer,
#     expectations=expectations,
#     overlaps=overlaps
# )
inv_obj_fit = inverse_objective_function(ising_model, N, p, mixer, None, None)
fitres = inv_obj_fit(np.hstack([gamma_fitted, beta_fitted]))
print("Fitted Proxy Result (negative expectation):", fitres)

# Store the expectations (note: inverse_objective_function returns -expectation)
initial_expectations.append(-initres)
fitted_expectations.append(-fitres)

initial_approx_ratios.append(maxcut_approx_ratio(graph, -initres))
fitted_approx_ratios.append(maxcut_approx_ratio(graph, -fitres))



# Print the mean expectations
print("Initial Proxy Mean Expectation:", np.mean(initial_expectations))
print("Fitted Proxy Mean Expectation :", np.mean(fitted_expectations))

print("Initial Proxy Mean Approx Ratio:", np.mean(initial_approx_ratios))
print("Fitted Proxy Mean Approx Ratio :", np.mean(fitted_approx_ratios))

print(f"Initial Proxy MSE: {initial_mse}")
print(f"Fitted  Proxy MSE: {fitted_mse}")
# %%
