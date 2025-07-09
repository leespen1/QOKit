#%%  imports 
#%% imports 
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import qokit.maxcut as mc
from grips.QAOA_simulator import get_expectation, get_simulator
from grips.QAOA_proxy_interface import QAOA_proxy, QAOA_proxy_expectation
import grips.triangle_proxy as tpr
import grips.paper_proxy as ppr
import grips.normal_proxy as npr
import os  
import grips.real_distribution as rd 
from grips.triangle_proxy import TriangleProxy

from juliacall import Main as jl
jl.seval('''
using Pkg
Pkg.activate(joinpath(@__DIR__, "../julia"))
Pkg.instantiate()
using JuliaQAOA
''')

# %%
num_edges = 10
edge_probability = 0.5
num_graphs = 10
graphs = [nx.erdos_renyi_graph(num_edges, edge_probability) for _ in range(num_graphs)]

realdist = rd.get_homogeneous_distribution(graphs)
print(realdist)
realdist.shape



# %%
import numpy as np
from scipy.optimize import minimize

def loss_function(params, realdist, num_constraints, num_qubits):
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
initial_params = [0.0, 0.0, 1.0, 1.0]  # reasonable starting point
num_constraints = num_edges  # Number of constraints
num_qubits = graphs[0].number_of_nodes()  # Number of qubits
bounds = [
    (0, None),   # h_tweak_sub >= 0
    (-10, 10),   # hc_tweak_add can be small positive/negative
    (0.1, 20),   # l_tweak_mul > 0
    (0.1, 20),   # r_tweak_mul > 0
]

result = minimize(
    loss_function,
    initial_params,
    args=(realdist, num_constraints, num_qubits),
    method='L-BFGS-B',
    bounds=bounds,
    options={'maxiter': 10000}
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
