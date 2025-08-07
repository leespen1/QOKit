import networkx as nx, qokit.maxcut as mc, numpy as np, typing
from numba import njit
from grips.QAOA_simulator import QAOA_run, get_simulator


def get_homogeneous_distribution(graph):
    """
    Given a graph, compute the real distribution n(x; d, c) (as a 3D array),
    then compute N(c'; d, c) by averaging n(x; d, c) over each x with cost c'.
    
    Given a LIST of graphs, uses all of them and averages their homogeneous dists 
    """
    if type(graph) is nx.Graph:
        costs = get_costs(graph)
        num_edges = graph.number_of_edges()
        num_vertices = graph.number_of_nodes()
        real_distribution = get_real_distribution_from_costs(costs, num_edges, num_vertices)[0]
        return get_homogeneous_distribution_from_costs(costs, real_distribution)
    elif type(graph) is list:
        costs = get_costs(graph)
        num_edges = graph[0].number_of_edges()
        num_vertices = graph[0].number_of_nodes()
        max_cost = max(g.number_of_edges() for g in graph)  # Max edges across all graphs
        #max_cost = np.floor(num_vertices/2)*np.ceil(num_vertices/2)  #alternate bound independent of edge number
        real_distributions = get_real_distribution_from_costs(costs, num_edges, num_vertices, max_possible_cost=max_cost)
        
        dist = None
        for i, realdist in enumerate(real_distributions):
            if dist is None:
                dist = get_homogeneous_distribution_from_costs(costs[i], realdist, max_possible_cost=max_cost)
            else:
                dist += get_homogeneous_distribution_from_costs(costs[i], realdist, max_possible_cost=max_cost)
        dist /= len(real_distributions)  # Average them
        return dist



def get_real_distribution(graph):
    """
    Given a graph, compute the real distribution n(x; d, c) from the paper.

    If a LIST of graphs is passed in, average the resulting distributions 
    for each graph in the list. 

    n(x; d, c) gives the number of bitstrings that have cost c and are Hamming
    distance c from bitstring x.

    The distribution is implemented as an 3D array:
    - Index 1 chooses the bitstring, e.g. 000, 001, 010, 011, etc
    - Index 2 chooses the Hamming distance, e.g. 0, 1, 2, ..., num_qubits
    - Index 3 chooses the cost, e.g. 0, 1, 2, ... 

    This file is for taking a particular graph, and obtaining the real distribution
    n(x; d, c) from the paper. That is, the number of bitstrings with cost c that
    are Hamming distance d from the bitstring x.
    """

    if type(graph) is nx.Graph:
        costs = get_costs(graph)
        num_edges = graph.number_of_edges()
        num_vertices = graph.number_of_nodes()
        return get_real_distribution_from_costs(costs, num_edges, num_vertices)[0]

    elif type(graph) is list:
        costs = get_costs(graph)
        num_edges = graph.number_of_edges()
        num_vertices = graph.number_of_nodes()
        max_cost = np.floor(num_vertices/2)*np.ceil(num_vertices/2)
        return get_real_distribution_from_costs(costs, num_edges, num_vertices, max_possible_cost=max_cost)


def get_costs(graph):
    """
    Given a graph, get the cost associated with each bitstring (for the maxcut
    problem), as a 1D array.
    If a LIST of graphs is passed as graph, return a 2D array where each 
    row is the costs for. 
    """

    if type(graph) is nx.Graph:
        num_vertices = graph.number_of_nodes()

        # Get the cost associated with each bitstring, as a 1D array
        ising_model = mc.get_maxcut_terms(graph)
        sim = get_simulator(num_vertices, ising_model)
        costs = sim.get_cost_diagonal()
        return costs
    
    elif type(graph) is list:
        num_graphs = len(graph)
        # If the graph is a list of graphs, compute the costs for each graph
        costs = np.zeros((num_graphs, 2 ** graph[0].number_of_nodes()))
        for i, g in enumerate(graph):
            costs[i] = get_costs(g) #SO CLEAN SELF-CALLING <3 (though we can make this faster actually)
        return costs
    
    else:
        raise TypeError(f"Argument 'graph' needs to be a networkx.Graph or a list of networkx.Graph, got {type(graph)} instead.")





@njit
def get_real_distribution_from_costs(costs, num_edges, num_vertices, max_possible_cost=None):
    """
    Given a the set of costs associated with each bitstring for a graph
    partition, the number of edges in the graph, and the number of vertices in
    the graph, compute the real distribution n(x; d, c) from the paper.

    Given a 2D np array of costs, returns a LIST of distritbutions. 
    n(x; d, c) gives the number of bitstrings that have cost c and are Hamming
    distance c from bitstring x.

    The distribution is implemented as an 3D array:
    - Index 1 chooses the bitstring, e.g. 000, 001, 010, 011, etc
    - Index 2 chooses the Hamming distance, e.g. 0, 1, 2, ..., num_qubits
    - Index 3 chooses the cost, e.g. 0, 1, 2, ... 

    This file is for taking a particular graph, and obtaining the real distribution
    n(x; d, c) from the paper. That is, the number of bitstrings with cost c that
    are Hamming distance d from the bitstring x.

    returns: numpy array n_distribution of shape (num_bitstrings, num_distances, num_costs), 
    with n_distribution[x, d, c] giving the number of bitstrings with cost c that are
    Hamming distance d from the bitstring x."""

    num_bitstrings = 2 ** num_vertices
    num_distances = 1 + num_vertices

    if max_possible_cost is None:
        num_costs = 1 + num_edges
        assert len(costs.shape) == 1
        assert len(costs) == num_bitstrings

        n_distribution = np.zeros((num_bitstrings, num_distances, num_costs), dtype = np.int32)

        for x in range(num_bitstrings): # 0:num_bitstrings-1
            for y in range(num_bitstrings):
                d = hamming_distance(x, y)
                cost_y = costs[y].astype(np.int32)  # Convert from float to int
                n_distribution[x, d, cost_y] += 1

        return [n_distribution]
    
    else:
        #this is for the case where we average multiple distributions
        #these have costs from 0 to max_possible_cost, and now costs is 
        # an array with shape(number of graphs, number of poss. costs)
        # costs has shape (number_of_graphs, num_bitstrings)
        assert len(costs.shape) == 2
        max_possible_cost = int(max_possible_cost) #should be max edge count of whole collection of graphs 
        num_graphs = int(costs.shape[0])
        num_costs = int(max_possible_cost + 1)
        n_distribution = np.zeros((num_bitstrings, num_distances, num_costs), dtype=np.int32)

        # Compute all pairwise Hamming distances once here, much faster 
        bitstrings = np.arange(num_bitstrings, dtype=np.uint32)
        xor = np.bitwise_xor(bitstrings[:, None], bitstrings[None, :])

        # Numba-compatible popcount implementation
        def popcount(x):
            count = 0
            while x:
                count += x & 1
                x >>= 1
            return count

        # Compute popcounts for all pairs
        popcounts = np.zeros_like(xor, dtype=np.int32)
        for i in range(xor.shape[0]):
            for j in range(xor.shape[1]):
                popcounts[i, j] = popcount(xor[i, j])

        # Precompute indexing arrays X and D (same for all graphs)
        X = np.repeat(np.arange(num_bitstrings), num_bitstrings)
        D = popcounts.ravel()

       
        n_distributions = [] #to return
        for i in range(num_graphs):
            n_dist_i = np.zeros((num_bitstrings, num_distances, num_costs), dtype=np.int32)

            costs_int = costs[i].astype(np.int32)
            C = np.empty(num_bitstrings * num_bitstrings, dtype=np.int32)
            for k in range(num_bitstrings):
                for l in range(num_bitstrings):
                    C[k * num_bitstrings + l] = costs_int[l]

            for idx in range(X.shape[0]):
                n_dist_i[X[idx], D[idx], C[idx]] += 1

            n_distributions.append(n_dist_i)

        return n_distributions



@njit
def get_homogeneous_distribution_from_costs(costs, real_distribution, max_possible_cost = None):
    """
    Given the costs associated with each bitstring, and the real distribution
    n(x; d, c) (as a 3D array), compute N(c'; d, c) by averaging n(x; d, c)
    over each x with cost c'.

    max_possible_cost: should be the number of edges, for the UNWEIGHTED maxcut problem.
    This is so that multiple homogeneous distributions can be averaged. 
    """

    if max_possible_cost is None:
        #this is spencer's version 
        assert len(real_distribution.shape) == 3
        num_bitstrings = real_distribution.shape[0]
        num_distances = real_distribution.shape[1]
        num_costs = real_distribution.shape[2]

        homogeneous_distribution = np.zeros((num_costs, num_distances, num_costs))
        num_cost_occurences = np.zeros(num_costs)

        # Sum over bitstrings with cost c'
        for bitstring_i in range(num_bitstrings):
            bitstring_cost = int(costs[bitstring_i])
            num_cost_occurences[bitstring_cost] += 1

            homogeneous_distribution[bitstring_cost, :, :] += real_distribution[bitstring_i, :, :]

    # Take average over number of bitstrings with cost c'
    for bitstring_cost in costs:
        bitstring_cost_int = int(bitstring_cost)
        homogeneous_distribution[bitstring_cost_int, :, :] /= num_cost_occurences[bitstring_cost_int]

    return homogeneous_distribution



@njit
def hamming_distance(bitstring1: int, bitstring2: int):
    """
    Compute the Hamming distance between two bitstrings, given as integers.

    The Hamming distance between bitstrings x and y is the number of bits which
    are different. E.g:
    - d(000, 111) = 3
    - d(100, 111) = 2
    - d(001, 111) = 2
    - d(101, 111) = 1
    - d(111, 111) = 0
    - d(010, 101) = 3
    """
    xor = bitstring1 ^ bitstring2
    d = bitcount(xor)
    # Bits of xor are 1 if-and-only-if the bits of x and y differ.
    # Therefore, xor.bit_count() is the number of bits that differ. 
    return d

@njit
def bitcount(x):
    count = 0
    while x:
        count += x & 1
        x >>= 1
    return count

