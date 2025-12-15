

import sys
import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import qokit.maxcut as mc

# Add grips to path
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..')))

from grips.proxies import PaperProxy
from grips.QAOA_proxy_interface import QAOA_proxy_optimize_gamma_beta, estimate_P_for_graphs
from grips import inverse_objective_function, maxcut_approx_ratio




# Parameters
GRAPH_TYPE = 'barabasi_albert' # 'erdos_renyi', 'barabasi_albert'
NUM_NODES = 12
NUM_GRAPHS = 1
# Graph specific params
BA_M = 2 # Edges to attach from a new node to existing nodes
ER_P = 0.5

# Sampling Params
SAMPLING_METHOD = 'moment_matching'
DIST_TYPE = 'edgeworth' # 'gaussian', 'beta', 'edgeworth'
NUM_SAMPLES = 1000

# QAOA Params
P_LAYERS = 1
OPTIMIZER_METHOD = 'COBYLA'
MAX_ITER = 100





# Generate Graphs
graphs = []
for _ in range(NUM_GRAPHS):
    if GRAPH_TYPE == 'barabasi_albert':
        g = nx.barabasi_albert_graph(NUM_NODES, BA_M)
    elif GRAPH_TYPE == 'erdos_renyi':
        g = nx.erdos_renyi_graph(NUM_NODES, ER_P)
    else:
        raise ValueError("Unknown graph type")
    graphs.append(g)

print(f"Generated {len(graphs)} {GRAPH_TYPE} graphs with {NUM_NODES} nodes.")





# Estimate Cost Distribution
print("Estimating cost distribution from samples...")
sampled_P = estimate_P_for_graphs(
    graphs, 
    method=SAMPLING_METHOD, 
    num_samples=NUM_SAMPLES, 
    dist_type=DIST_TYPE
)
print("Done.")





# Initialize Proxies
num_edges_list = [g.number_of_edges() for g in graphs]
avg_edges = int(np.mean(num_edges_list))
print(f"Average number of edges: {avg_edges}")

# Default Paper Proxy (uses Binomial distribution internally)
# PaperProxy needs num_constraints, num_qubits, prob_edge.
# For BA graphs, prob_edge isn't well defined like ER. 
# We can estimate effective prob_edge = 2*E / (N*(N-1))
prob_edge_eff = 2 * avg_edges / (NUM_NODES * (NUM_NODES - 1))
print(f"Effective edge probability: {prob_edge_eff:.4f}")

default_proxy = PaperProxy(
    num_constraints=avg_edges,
    num_qubits=NUM_NODES,
    prob_edge=prob_edge_eff
)

# Sampled Paper Proxy
sampled_proxy = PaperProxy(
    num_constraints=avg_edges,
    num_qubits=NUM_NODES,
    prob_edge=prob_edge_eff,
    cost_distribution=sampled_P
)




# Optimize QAOA Parameters
print("Optimizing parameters for Default Proxy...")
init_gamma = np.array([0.1] * P_LAYERS)
init_beta = np.array([0.1] * P_LAYERS)

res_default = QAOA_proxy_optimize_gamma_beta(
    default_proxy, 
    init_gamma, 
    init_beta, 
    optimizer_method=OPTIMIZER_METHOD,
    optimizer_options={'maxiter': MAX_ITER}
)
gamma_default = res_default['gamma']
beta_default = res_default['beta']
print(f"Default Proxy Params: Gamma={gamma_default}, Beta={beta_default}")

print("\nOptimizing parameters for Sampled Proxy...")
res_sampled = QAOA_proxy_optimize_gamma_beta(
    sampled_proxy, 
    init_gamma, 
    init_beta, 
    optimizer_method=OPTIMIZER_METHOD,
    optimizer_options={'maxiter': MAX_ITER}
)
gamma_sampled = res_sampled['gamma']
beta_sampled = res_sampled['beta']
print(f"Sampled Proxy Params: Gamma={gamma_sampled}, Beta={beta_sampled}")





# Evaluate on Real Graphs
print("\nEvaluating on real graphs...")

default_ratios = []
sampled_ratios = []

for g in graphs:
    ising_model = mc.get_maxcut_terms(g)
    N = g.number_of_nodes()
    
    # Function to get expectation value (negative of cost)
    # Using python simulator because ofLLVM/Numba issues with the default causing crash...?
    obj_func = inverse_objective_function(ising_model, N, P_LAYERS, "x", None, None, simulator_name="python")
    
    # Evaluate Default Params
    neg_exp_default = obj_func(np.hstack([gamma_default, beta_default]))
    exp_default = -neg_exp_default
    ratio_default = maxcut_approx_ratio(g, exp_default)
    default_ratios.append(ratio_default)
    
    # Evaluate Sampled Params
    neg_exp_sampled = obj_func(np.hstack([gamma_sampled, beta_sampled]))
    exp_sampled = -neg_exp_sampled
    ratio_sampled = maxcut_approx_ratio(g, exp_sampled)
    sampled_ratios.append(ratio_sampled)

print(f"Mean Approx Ratio (Default Proxy): {np.mean(default_ratios):.4f}")
print(f"Mean Approx Ratio (Sampled Proxy): {np.mean(sampled_ratios):.4f}")



# Plotting
plt.figure(figsize=(10, 6))
plt.boxplot([default_ratios, sampled_ratios], labels=['Default Proxy', 'Sampled Proxy'])
plt.title(f'QAOA Performance Comparison (p={P_LAYERS})\n{GRAPH_TYPE} graphs, N={NUM_NODES}')
plt.ylabel('Approximation Ratio')
plt.grid(True, axis='y')
plt.show()