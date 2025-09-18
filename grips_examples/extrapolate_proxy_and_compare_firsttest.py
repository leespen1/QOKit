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
    get_homogeneous_distribution,get_homogeneous_distribution_from_proxy,
    TriangleProxy, NormalProxy, PaperProxy,
    inverse_objective_function, get_expectation,  
    maxcut, maxcut_approx_ratio, spsa_for_scipy,
    plot_distribution_lines_all, fit_proxy_to_real, 
    pad_and_stack, jl
)

from scipy.optimize import minimize
from scipy.optimize import dual_annealing
print("Finished importing python packages/functions!")


#%% Params to alternate against initial

'''
probably should turn distfit_multiple_graphs into a 
function so we can call it here to get the fitted parameters. 
These were just obtained manually to test: 

With 5 nodes, 100 graphs, edge probability 0.5, 
Fitted parameters were [1.0, 1.11, 2.31, 0.1]

With 6 nodes, 100 graphs, edge probability 0.5, 
Fitted parameters were [1.33, 1.11, 2.31, 0.1]

With 7 nodes, 100 graphs, edge probability 0.5, 
Fitted parameters were [4.57, 2.86, 2.94, 0.1]

With 8 nodes, 100 graphs, edge probability 0.5, 
Fitted parameters were [9.33, 1.11, 2.31, 0.1]

With 9 nodes, 100 graphs, edge probability 0.5, 
Fitted parameters were [21.88, 2.5, 2.59, 0.1]

With 10 nodes, 100 graphs, edge probability 0.5, 
Fitted parameters were [35.0, -2.86, 1.05, 0.1]
 '''
data = {5:[1.0, 1.11, 2.31, 0.1], 
        6:[1.33, 1.11, 2.31, 0.1],
        7:[4.57, 2.86, 2.94, 0.1],
        8:[9.33, 1.11, 2.31, 0.1],
        9:[21.88, 2.5, 2.59, 0.1],
        10:[35.0, -2.86, 1.05, 0.1]
        }

# Fit, plot, and predict parameters
import numpy as np
import matplotlib.pyplot as plt

# Extract data for fitting
nodes = np.array(list(data.keys()))
params = np.array(list(data.values()))
first_param = params[:, 0]

# Perform quadratic fit for the first parameter
# f(x) = c*x**2 + b*x + a
coeffs = np.polyfit(nodes, first_param, 2)
poly = np.poly1d(coeffs)

# Calculate mean for other parameters
mean_other_params = np.mean(params[:, 1:], axis=0)

# Predict parameters for 11 nodes
nodes_to_predict = 11
pred_first_param = poly(nodes_to_predict)
params_11_nodes = [pred_first_param] + list(mean_other_params)

print(f"Quadratic fit for h_tweak_sub: f(x) = {coeffs[0]:.2f}x^2 + {coeffs[1]:.2f}x + {coeffs[2]:.2f}")
print(f"Mean of other parameters: {mean_other_params}")
print(f"Predicted parameters for {nodes_to_predict} nodes: {params_11_nodes}")

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(nodes, first_param, 'o', label='Actual Data')
x_fit = np.linspace(min(nodes), nodes_to_predict, 100)
y_fit = poly(x_fit)
plt.plot(x_fit, y_fit, '-', label='Quadratic Fit')
plt.xlabel("Number of Nodes")
plt.ylabel("First Parameter Value")
plt.title("Quadratic Fit of First Parameter vs. Number of Nodes")
plt.legend()
plt.grid(True)
plt.show()

#%%
# %% graphs
print("\nSetting up a graphs, getting dist...")
num_nodes = 6
num_graphs = 10
edge_probability = 0.4
graphs = []

for i in range(num_graphs):
    graph = nx.erdos_renyi_graph(num_nodes, edge_probability) 
    graphs.append(graph)



initial_params = [0, 0, 1, 1]  
alternate_params = [0.5, 0, 2, 2]  
max_num_edges = max(g.number_of_edges() for g in graphs)
num_constraints = max_num_edges # Number of constraints -- +1 here caused error previously!
num_qubits = num_nodes  # Number of qubits

# proxy = TriangleProxy(num_constraints, num_qubits, *initial_params)
# initial_triangle_homodist = get_homogeneous_distribution_from_proxy(proxy)
# plot_distribution_lines_all(initial_triangle_homodist, f"Initial Triangle Proxy Distribution {initial_params}")

# %% Define mean-squared error loss function for fitting triangle proxy to statistical homogeneous distribution
print("\nDefining triangle proxy loss function ... ")

def mse_dist_loss(params, homodist, num_constraints=0):
    h_tweak_sub, hc_tweak_add, l_tweak_mul, r_tweak_mul = params
    num_constraints = max(homodist.shape[0] - 1, num_constraints)
    num_qubits = homodist.shape[1] - 1

    proxy = TriangleProxy(
        num_constraints=num_constraints,
        num_qubits=num_qubits,
        h_tweak_sub=h_tweak_sub,
        hc_tweak_add=hc_tweak_add,
        l_tweak_mul=l_tweak_mul,
        r_tweak_mul=r_tweak_mul
    )
    return grips.distribution_mean_squared_error(proxy, homodist)

print("Finished defining triangle proxy loss function!")

#%%

print("alternate parameters:", alternate_params)
alternate_proxy = TriangleProxy(
    num_constraints=num_constraints,
    num_qubits=num_qubits,
    h_tweak_sub=alternate_params[0],
    hc_tweak_add=alternate_params[1],
    l_tweak_mul=alternate_params[2],
    r_tweak_mul=alternate_params[3]
)



#%% Run QAOA with alternate proxy
print("\nRunning QAOA with alternate proxy (but not tuned gamma/beta)...")
gammas = np.linspace(0, np.pi, 10)  # Gamma values for QAOA
betas = np.linspace(0, np.pi, 10)  # Beta values for QAOA
alternate_triangle_results = QAOA_proxy(alternate_proxy, gammas, betas)
print("alternate Proxy Results:", alternate_triangle_results)

final_amplitudes = alternate_triangle_results[-1]
expectation = QAOA_proxy_expectation(alternate_proxy, final_amplitudes)
print("Expectation value:", expectation)
print("Finished running QAOA with alternate proxy!")

# %% Run QAOA with initial, unalternate proxy
print("\nRunning QAOA with initial, unalternate proxy ...")
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

# # alternate expectations
# print(f"Initial Expectation: {initial_expectation}")
# print(f"alternate Expectation: {expectation}")
print("Finished running QAOA with initial, unalternate proxy!")


# Comparing QAOA results for initial versus alternate proxy
'''
-defining initial gammas and betas
-finding best gamma and beta according to initial versus alternate proxy
-using these to run QAOA and comparing expectations
'''

gamma_0 = np.array([0.1])
beta_0 = np.array([0.1])
init_result = QAOA_proxy_optimize_gamma_beta(initial_proxy, gamma_0, beta_0, optimizer_method = 'Nelder-Mead', optimizer_options={'maxiter': 1})
gamma_init = init_result["gamma"]
beta_init = init_result["beta"]



alternate_result = QAOA_proxy_optimize_gamma_beta(alternate_proxy, gamma_0, beta_0, optimizer_method = 'Nelder-Mead', 
                                               optimizer_options = {'maxiter': 10000, 'epsilon': 0.00001})
gamma_alternate = alternate_result["gamma"]
beta_alternate = alternate_result["beta"]

graph = nx.erdos_renyi_graph(num_nodes, edge_probability)


#%% get QAOA expectations of inverse objective function for initial and alternate proxies
initial_expectations = []
alternate_expectations = []
initial_approx_ratios = []
alternate_approx_ratios = []

for graph in graphs:
    ising_model = mc.get_maxcut_terms(graph)
    N = graph.number_of_nodes()
    sim = get_simulator(N, ising_model)
    p = 1
    mixer = "x"
    expectations = []
    overlaps = []

    '''
    we want the QAOA expectations of the alternate gammas and bets from init proxy versus alternate proxy, 
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
    # print("Initial Proxy Result (negative expectation):", initres)

    # For alternate parameters
    expectations = []
    overlaps = []
    # fitres = QAOA_run(
    #     ising_model=ising_model,
    #     N=N,
    #     p=p,
    #     init_gamma=gamma_alternate,
    #     init_beta=beta_alternate,
    #     optimizer_method='COBYLA',
    #     optimizer_options=None,
    #     mixer=mixer,
    #     expectations=expectations,
    #     overlaps=overlaps
    # )
    inv_obj_fit = inverse_objective_function(ising_model, N, p, mixer, None, None)
    fitres = inv_obj_fit(np.hstack([gamma_alternate, beta_alternate]))
    # print("alternate Proxy Result (negative expectation):", fitres)

    # Store the expectations (note: inverse_objective_function returns -expectation)
    initial_expectations.append(-initres)
    alternate_expectations.append(-fitres)

    initial_approx_ratios.append(maxcut_approx_ratio(graph, -initres))
    alternate_approx_ratios.append(maxcut_approx_ratio(graph, -fitres))


# Print the mean expectations
print("\n\nInitial Proxy Mean Expectation:", np.mean(initial_expectations))
print("alternate Proxy Mean Expectation :", np.mean(alternate_expectations))

print("\n\nInitial Proxy Mean Approx Ratio:", np.mean(initial_approx_ratios))
print("alternate Proxy Mean Approx Ratio :", np.mean(alternate_approx_ratios))

