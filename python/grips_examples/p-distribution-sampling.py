import sys
import os
import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import networkx as nx
from grips.QAOA_proxy_interface import estimate_P_for_graphs


# --- Example Usage of cost distribution estimation
graphs = [nx.erdos_renyi_graph(50, 0.5) for _ in range(3)]

# The 'beta' distribution is generally superior for this problem class
P_beta = estimate_P_for_graphs(graphs, num_samples=10000, method="moment_matching", dist_type="edgeworth")

# Plotting
mu, sigma, _, _ = P_beta.params
x = np.linspace(mu - 4 * sigma, mu + 4 * sigma, 1000)
y = [P_beta(val) for val in x]

plt.figure(figsize=(10, 6))
plt.plot(x, y, label="Edgeworth Approximation")
plt.title("Cost Distribution Approximation")
plt.xlabel("Cut Size")
plt.ylabel("Probability Density")
plt.legend()
plt.grid(True)
plt.savefig("cost_distribution.png")
print("Plot saved to cost_distribution.png")
