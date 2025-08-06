################################################################################
#
# test_real_distribution_2_node_graph.py
#
# This script tests that the function `get_real_distribution` produces the
# correct real distribution n(c';d,c) for the simple, 2-node graph below:
#
#      [0]---[1]
#
# Note that the nodes MUST be numbered 0 and 1. If they are numbered as 1 and 2,
# then the graph will be interpreted as
#
#      [0]    [1]---[2]
#
# In that case, the correspondence between the bitstrings and the graph
# partitions is different. For example, the bitstring 0b001 places the "0" node
# in its own partition. This results in a cost of 1 for the top graph, and 0
# for the bottom graph.
#
################################################################################
from grips.real_distribution import get_real_distribution, get_homogeneous_distribution
import networkx as nx, numpy as np

# Create a graph that consists of two connected nodes
G = nx.Graph()
G.add_edge(0,1)

n_dist = get_real_distribution(G)

print("Computed n_dist:")
print(n_dist.astype(int))

correct_n_dist = np.array([
  [[1, 0],
   [0, 2],
   [1, 0]],

  [[0,1],
   [2,0],
   [0,1]],

  [[0,1],
   [2,0],
   [0,1]],

  [[1,0],
   [0,2],
   [1,0]]
])

assert np.array_equal(n_dist.astype(int), correct_n_dist)
print("Real distribution n_dist was computed correctly for the two-node example.")

N_dist = get_homogeneous_distribution(G)

print("Computed N_dist:")
print(N_dist.astype(int))

correct_N_dist = np.array([
  [[1, 0],
   [0, 2],
   [1, 0]],

  [[0,1],
   [2,0],
   [0,1]],
])

assert np.array_equal(N_dist.astype(int), correct_N_dist)
print("Homogenous distribution N_dist was computed correctly for the two-node example.")
