import networkx as nx
# import dwave_networkx as dnx
# import itertools


import networkx as nx
from typing import Tuple, Set, Any, Optional

def maxcut(G: nx.Graph) -> Tuple[int, Optional[Set[Any]]]:
    """
    Solve the Max-Cut by brute force.
    Args:
        G: A networkx Graph object.

    Returns:
        - objective_value (int): The value of the maximum cut.
        - best_partition (set): One of the two sets in the optimal partition
                                achieving the max cut.
    """
    nodes = list(G.nodes)
    n = len(nodes)

    assert n != 0, "Your input graph is empty!"

    max_cut_value = 0
    best_S = set() 

    # Brute force! (note: don't need to check symmetric partitions)
    for i in range(2**(n-1)):
        S = {nodes[0]}  # Fix the first node to be in S.
        
        # Iterate through the remaining n-1 nodes.
        for j in range(n - 1):
            # Check the j-th bit of i to decide which partition to place node[j+1] in.
            if (i >> j) & 1:
                S.add(nodes[j + 1])
        
        #get the cut value
        cut_value = nx.cut_size(G, S)

        if cut_value > max_cut_value:
            max_cut_value = cut_value
            best_S = S

    return max_cut_value, best_S

def maxcut_approx_ratio(G, approx_objective_value):
    """
    Computes the approximation ratio of an approximate MaxCut solution.
    Args:
        G: networkx graph
        approx_objective_value: objective value of approximate solution
    Returns:
        ratio: approx_objective_value / optimal_objective_value
    """
    optimal_objective_value, _ = maxcut(G)
    assert optimal_objective_value != 0,\
        'Error: Optimal objective value is zero!!'

    return approx_objective_value / optimal_objective_value