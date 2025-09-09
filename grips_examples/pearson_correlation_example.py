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
    pad_and_stack, get_pearson_correlation_coefficients, pad_to_match
)

def n_choose_2(n):
    return n * (n - 1) // 2

# 10 vertices, edge probability 1/3, same as Fig 3 of paper
num_nodes = 10
num_graphs = 100
edge_probability = 1/3
graphs = []
homodists = []
import random
for i in range(num_graphs):
    graph = nx.erdos_renyi_graph(num_nodes, edge_probability)
    #approx ratios later break if we have graphs with no edges
    #add a random edge if the graph doesn't have one
    if graph.number_of_edges() == 0:
        nodes = list(graph.nodes)
        if len(nodes) >= 2:
            u, v = random.sample(nodes, 2)
            graph.add_edge(u, v)
    homodist = get_homogeneous_distribution(graph)
    graphs.append(graph)
    homodists.append(homodist)

#Q: is the scaling correct here with np.mean()?
homodists = pad_and_stack(homodists)
analytic_homodist = np.mean(homodists, axis=0)

expected_number_of_edges = int(edge_probability * n_choose_2(num_nodes))
paper_proxy = PaperProxy(expected_number_of_edges, num_nodes, edge_probability)
triangle_proxy = TriangleProxy(expected_number_of_edges, num_nodes)
paper_homodist = get_homogeneous_distribution_from_proxy(paper_proxy)
triangle_homodist = get_homogeneous_distribution_from_proxy(triangle_proxy)

# pad shapes to match
analytic_homodist, paper_homodist = pad_to_match(analytic_homodist, paper_homodist)
analytic_homodist, triangle_homodist = pad_to_match(analytic_homodist, triangle_homodist)

paper_pearson_coefficients = get_pearson_correlation_coefficients(analytic_homodist, paper_homodist)
triangle_pearson_coefficients = get_pearson_correlation_coefficients(analytic_homodist, triangle_homodist)


plt.plot(range(analytic_homodist.shape[0]), paper_pearson_coefficients, label="Paper Proxy", marker='o')
plt.plot(range(analytic_homodist.shape[0]), triangle_pearson_coefficients, label="Triangle Proxy", marker='o')
plt.ylim(0, 1.1)
plt.grid(True, which="both")
plt.xlabel("Cost (c') of Target Bitstring")
plt.ylabel("Pearson Correlation Coefficient")
plt.show()

